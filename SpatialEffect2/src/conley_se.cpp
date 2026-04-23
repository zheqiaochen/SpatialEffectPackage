#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Conley spatial HAC standard errors
//
// Computes the meat of the sandwich variance estimator with spatial kernel weighting.
// When cutoff == 0, returns EHW (heteroskedasticity-robust) meat instead.
//
// @param res       Residual vector (N x 1)
// @param W_meat    Design matrix for meat (N x K), possibly pre-weighted
// @param dist      Pre-computed distance matrix (N x N)
// @param XXinv     (X'X)^{-1} matrix (K x K)
// @param cutoff    Spatial bandwidth for kernel weighting (0 = EHW)
// @param kernel    Kernel type: 1=uniform, 2=triangular, 3=Epanechnikov
// @param trim      Whether to trim kernel matrix eigenvalues for PSD
// @param if_edof   Whether to compute effective degrees of freedom
// @param n_threads Number of OpenMP threads
//
// Returns: list(VCE_meat, Dist_kernel, mu, v)
// [[Rcpp::export]]
Rcpp::List ConleySE(const arma::vec& res, const arma::mat& W_meat,
                    const arma::mat& dist, const arma::mat& XXinv,
                    double cutoff, int kernel, int trim, int if_edof,
                    int n_threads = 1) {
  int N = res.n_elem;
  int K = W_meat.n_cols;
  double mu = 1.0;
  double v = 2.0;

  // Build kernel matrix
  arma::mat dist_kernel(N, N, arma::fill::zeros);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  if (cutoff > 0) {
#pragma omp parallel for schedule(static) if(N > 100)
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        double u = dist(i, j) / cutoff;
        if (u <= 1.0) {
          if (kernel == 3) {
            dist_kernel(i, j) = 0.75 * (1.0 - u * u);
          } else if (kernel == 1) {
            dist_kernel(i, j) = 1.0;
          } else {
            dist_kernel(i, j) = 1.0 - u;
          }
        }
      }
    }
  }

  // Eigenvalue trimming
  if (trim == 1 && cutoff > 0) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, dist_kernel);
    eigval.elem(arma::find(eigval < 0)).zeros();
    dist_kernel = eigvec * arma::diagmat(eigval) * eigvec.t();
  }

  // Effective degrees of freedom (Welch-Satterthwaite)
  if (if_edof == 1) {
    arma::mat M_mat = arma::eye(N, N) - W_meat * XXinv * W_meat.t();
    arma::vec w_coef = arma::zeros(K);
    w_coef(1) = 1.0;
    arma::vec z_vec = W_meat * XXinv * w_coef;

    arma::mat B(N, N, arma::fill::zeros);
#pragma omp parallel for schedule(dynamic) if(N > 100)
    for (int i = 0; i < N; i++) {
      B.col(i) = (z_vec % dist_kernel.col(i)) / XXinv(1, 1) * z_vec(i)
        - W_meat * XXinv * (W_meat.t() * (z_vec % dist_kernel.col(i))) / XXinv(1, 1) * z_vec(i);
    }
    mu = arma::trace(B);
    v = 2.0 * arma::trace(B * B);
  }

  // Compute meat matrix
  arma::mat meat(K, K, arma::fill::zeros);
  arma::mat meat_EHW(K, K, arma::fill::zeros);

  // Use parallelized reduction
  if (cutoff > 0) {
#pragma omp parallel if(N > 100)
    {
      arma::mat meat_local(K, K, arma::fill::zeros);
#pragma omp for schedule(static)
      for (int i = 0; i < N; i++) {
        arma::rowvec ri_Xi = res(i) * W_meat.row(i);
        arma::vec kernel_res = res % dist_kernel.col(i);
        meat_local += ri_Xi.t() * (kernel_res.t() * W_meat);
      }
#pragma omp critical
      meat += meat_local;
    }
  }

  // EHW (always computed when cutoff == 0)
#pragma omp parallel if(N > 100)
  {
    arma::mat ehw_local(K, K, arma::fill::zeros);
#pragma omp for schedule(static)
    for (int i = 0; i < N; i++) {
      ehw_local += res(i) * res(i) * W_meat.row(i).t() * W_meat.row(i);
    }
#pragma omp critical
    meat_EHW += ehw_local;
  }

  if (cutoff == 0) {
    return Rcpp::List::create(
      Rcpp::Named("VCE_meat") = meat_EHW,
      Rcpp::Named("Dist_kernel") = dist_kernel,
      Rcpp::Named("mu") = mu,
      Rcpp::Named("v") = v
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("VCE_meat") = meat,
      Rcpp::Named("Dist_kernel") = dist_kernel,
      Rcpp::Named("mu") = mu,
      Rcpp::Named("v") = v
    );
  }
}

// Variant that takes a pre-computed kernel matrix
// [[Rcpp::export]]
Rcpp::List ConleySE2(const arma::vec& res, const arma::mat& W_meat,
                     const arma::mat& dist_kernel, const arma::mat& XXinv,
                     int if_edof, int n_threads = 1) {
  int N = res.n_elem;
  int K = W_meat.n_cols;
  double mu = 1.0;
  double v = 2.0;

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  if (if_edof == 1) {
    arma::vec w_coef = arma::zeros(K);
    w_coef(1) = 1.0;
    arma::vec z_vec = W_meat * XXinv * w_coef;

    arma::mat B(N, N, arma::fill::zeros);
#pragma omp parallel for schedule(dynamic) if(N > 100)
    for (int i = 0; i < N; i++) {
      B.col(i) = (z_vec % dist_kernel.col(i)) / XXinv(1, 1) * z_vec(i)
        - W_meat * XXinv * (W_meat.t() * (z_vec % dist_kernel.col(i))) / XXinv(1, 1) * z_vec(i);
    }
    mu = arma::trace(B);
    v = 2.0 * arma::trace(B * B);
  }

  arma::mat meat(K, K, arma::fill::zeros);
#pragma omp parallel if(N > 100)
  {
    arma::mat meat_local(K, K, arma::fill::zeros);
#pragma omp for schedule(static)
    for (int i = 0; i < N; i++) {
      arma::rowvec ri_Xi = res(i) * W_meat.row(i);
      arma::vec kernel_res = res % dist_kernel.col(i);
      meat_local += ri_Xi.t() * (kernel_res.t() * W_meat);
    }
#pragma omp critical
    meat += meat_local;
  }

  return Rcpp::List::create(
    Rcpp::Named("VCE_meat") = meat,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("v") = v
  );
}

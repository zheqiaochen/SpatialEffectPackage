#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Compute kernel-weighted distance matrix
// kernel: 1 = uniform, 2 = triangular, 3 = Epanechnikov
// trim: 1 = eigenvalue trimming for positive semidefiniteness
// [[Rcpp::export]]
arma::mat KernelMatrix(const arma::mat& dist, double cutoff, int kernel = 1,
                       int trim = 0, int n_threads = 1) {
  int N = dist.n_rows;
  arma::mat K(N, N, arma::fill::zeros);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

#pragma omp parallel for schedule(static) if(N > 100)
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double u = dist(i, j) / cutoff;
      if (u <= 1.0) {
        if (kernel == 3) {
          K(i, j) = 0.75 * (1.0 - u * u);
        } else if (kernel == 1) {
          K(i, j) = 1.0;
        } else {
          K(i, j) = 1.0 - u;
        }
      }
    }
  }

  // Eigenvalue trimming for positive semidefiniteness
  if (trim == 1) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, K);
    eigval.elem(arma::find(eigval < 0)).zeros();
    K = eigvec * arma::diagmat(eigval) * eigvec.t();
  }

  return K;
}

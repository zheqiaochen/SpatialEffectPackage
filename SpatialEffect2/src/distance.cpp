#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Compute N x N Euclidean distance matrix (symmetric, parallelized)
// [[Rcpp::export]]
arma::mat DistMatEuclidean(const arma::vec& x, const arma::vec& y, int n_threads = 1) {
  int N = x.n_elem;
  arma::mat dist(N, N, arma::fill::zeros);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

#pragma omp parallel for schedule(dynamic) if(N > 100)
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      double dx = x(i) - x(j);
      double dy = y(i) - y(j);
      double d = std::sqrt(dx * dx + dy * dy);
      dist(i, j) = d;
      dist(j, i) = d;
    }
  }
  return dist;
}

// Compute Ny x Nz rectangular Euclidean distance matrix (parallelized)
// For distances between outcome locations and intervention nodes
// [[Rcpp::export]]
arma::mat DistMatEuclidean2(const arma::vec& x1, const arma::vec& y1,
                            const arma::vec& x2, const arma::vec& y2,
                            int n_threads = 1) {
  int Ny = x1.n_elem;
  int Nz = x2.n_elem;
  arma::mat dist(Ny, Nz);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

#pragma omp parallel for schedule(dynamic) if(Ny > 100)
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nz; j++) {
      double dx = x1(i) - x2(j);
      double dy = y1(i) - y2(j);
      dist(i, j) = std::sqrt(dx * dx + dy * dy);
    }
  }
  return dist;
}

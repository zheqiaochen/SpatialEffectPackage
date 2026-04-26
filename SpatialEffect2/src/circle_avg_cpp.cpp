#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Fast vectorized circle average computation for Euclidean point interventions
// on raster data. Instead of looping per node in R, this does everything in C++.
//
// @param raster_vals  Raster values as a flat vector (row-major)
// @param ncol_ras     Number of columns in the raster
// @param nrow_ras     Number of rows in the raster
// @param xmin         Raster extent x-minimum
// @param ymax         Raster extent y-maximum
// @param xres         Raster x-resolution
// @param yres         Raster y-resolution
// @param centers_x    X-coordinates of intervention nodes
// @param centers_y    Y-coordinates of intervention nodes
// @param offsets_x    X-offsets for sampling circle template
// @param offsets_y    Y-offsets for sampling circle template
// @param only_unique  Whether to deduplicate grid cells per node
// @param n_threads    Number of OpenMP threads
//
// Returns: list(Ybard, Ybard_sum, Ybard_len) each of length nz
// [[Rcpp::export]]
Rcpp::List CircleAvgRaster(const arma::vec& raster_vals,
                           int ncol_ras, int nrow_ras,
                           double xmin, double ymax,
                           double xres, double yres,
                           const arma::vec& centers_x,
                           const arma::vec& centers_y,
                           const arma::vec& offsets_x,
                           const arma::vec& offsets_y,
                           int only_unique,
                           int n_threads) {
  int nz = centers_x.n_elem;
  int n_offsets = offsets_x.n_elem;
  int ncells = ncol_ras * nrow_ras;

  arma::vec Ybard(nz);
  arma::vec Ybard_sum(nz, arma::fill::zeros);
  arma::ivec Ybard_len(nz, arma::fill::zeros);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

#pragma omp parallel for schedule(dynamic) if(nz > 10)
  for (int i = 0; i < nz; i++) {
    double cx = centers_x(i);
    double cy = centers_y(i);
    double sum_val = 0.0;
    int count = 0;

    // Track visited cells for deduplication
    std::vector<int> visited;
    if (only_unique == 1) {
      visited.reserve(n_offsets);
    }

    for (int k = 0; k < n_offsets; k++) {
      double px = cx + offsets_x(k);
      double py = cy + offsets_y(k);

      // Convert to raster cell index (row, col)
      int col = (int)std::floor((px - xmin) / xres);
      int row = (int)std::floor((ymax - py) / yres);

      if (col < 0 || col >= ncol_ras || row < 0 || row >= nrow_ras) continue;

      int cell_idx = row * ncol_ras + col;
      if (cell_idx < 0 || cell_idx >= ncells) continue;

      // Deduplication check
      if (only_unique == 1) {
        bool found = false;
        for (size_t v = 0; v < visited.size(); v++) {
          if (visited[v] == cell_idx) { found = true; break; }
        }
        if (found) continue;
        visited.push_back(cell_idx);
      }

      double val = raster_vals(cell_idx);
      if (!std::isnan(val)) {
        sum_val += val;
        count++;
      }
    }

    Ybard_sum(i) = sum_val;
    Ybard_len(i) = count;
    Ybard(i) = (count > 0) ? sum_val / count : R_NaN;
  }

  return Rcpp::List::create(
    Rcpp::Named("Ybard") = Ybard,
    Rcpp::Named("Ybard_sum") = Ybard_sum,
    Rcpp::Named("Ybard_len") = Ybard_len
  );
}

// Fast circle average for kriging predictions
// Evaluates kriging at sampling points and averages per node
// This version takes pre-computed sampling coordinates grouped by node
//
// @param predictions  Pre-computed kriging predictions at all sampling points
// @param group_ids    Node index (0-based) for each sampling point
// @param nz           Number of intervention nodes
//
// Returns: list(Ybard, Ybard_sum, Ybard_len)
// [[Rcpp::export]]
Rcpp::List CircleAvgFromPredictions(const arma::vec& predictions,
                                    const arma::ivec& group_ids,
                                    int nz) {
  arma::vec Ybard(nz);
  arma::vec Ybard_sum(nz, arma::fill::zeros);
  arma::ivec Ybard_len(nz, arma::fill::zeros);

  for (int k = 0; k < (int)predictions.n_elem; k++) {
    int g = group_ids(k);
    if (g >= 0 && g < nz && !std::isnan(predictions(k))) {
      Ybard_sum(g) += predictions(k);
      Ybard_len(g) += 1;
    }
  }

  for (int i = 0; i < nz; i++) {
    Ybard(i) = (Ybard_len(i) > 0) ? Ybard_sum(i) / Ybard_len(i) : R_NaN;
  }

  return Rcpp::List::create(
    Rcpp::Named("Ybard") = Ybard,
    Rcpp::Named("Ybard_sum") = Ybard_sum,
    Rcpp::Named("Ybard_len") = Ybard_len
  );
}

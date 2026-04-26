# Cross-validation for bandwidth selection

Two-stage grid search with K-fold (optionally spatial block)
cross-validation.

## Usage

``` r
CrossValidation(
  Sdata,
  outcome,
  treatment,
  dVec,
  grid = NULL,
  nfold = 5,
  block_cv = TRUE,
  parallel = TRUE,
  metric = "MSPE",
  kernel = "uni",
  bias_correction = FALSE,
  Zdata = NULL,
  x_coord_Z = NULL,
  y_coord_Z = NULL,
  n_threads = 1L
)
```

## Arguments

- Sdata:

  Stacked data frame

- outcome, treatment:

  Column name strings

- dVec:

  Distance vector

- grid:

  Bandwidth grid (auto-generated if NULL)

- nfold:

  Number of folds

- block_cv:

  Use spatial block CV

- parallel:

  Use future.apply parallelism

- metric:

  "MSPE" or "MAPE"

- kernel:

  Kernel type

- bias_correction:

  Use bias correction

- Zdata, x_coord_Z, y_coord_Z:

  Intervention data (required for block_cv)

- n_threads:

  C++ threads

## Value

List with CV.out, opt.bw, opt.bw.debias

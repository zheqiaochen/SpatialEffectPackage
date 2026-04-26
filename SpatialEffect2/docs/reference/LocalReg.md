# Local polynomial regression for smoothing AME estimates

Implements the local linear regression smoother with optional bias
correction via higher-order polynomial terms, following Hainmueller et
al. (2018).

## Usage

``` r
LocalReg(
  dVec,
  Sdata,
  bw,
  bw_debias = NULL,
  Zup = NULL,
  xevals = NULL,
  smooth.conley.se = 1,
  kernel = "uni",
  cutoff = 0,
  dist = NULL,
  dist.metric = "Euclidean",
  bias_correction = TRUE,
  Zdata = NULL,
  x_coord_Z = NULL,
  y_coord_Z = NULL,
  n_threads = 1L
)
```

## Arguments

- dVec:

  Distance values

- Sdata:

  Stacked data frame with outcome, treatment, dVec, pweight columns

- bw:

  Bandwidth

- bw_debias:

  Bandwidth for bias correction polynomial

- Zup:

  Treatment vector

- xevals:

  Evaluation points (defaults to dVec)

- smooth.conley.se:

  Compute Conley SE for smoothed estimates

- kernel:

  Kernel type string

- cutoff:

  Conley SE bandwidth

- dist:

  Pre-computed distance matrix

- dist.metric:

  Distance metric

- bias_correction:

  Use bias correction

- Zdata:

  Intervention data frame

- x_coord_Z, y_coord_Z:

  Coordinate column names

- n_threads:

  Number of threads

## Value

List with coefs, coefs_all, ses, coefs_debias

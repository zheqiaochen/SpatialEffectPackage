# Estimate the Average Marginalized Effect (AME) Curve

Estimates how treatment effects vary with distance in space for spatial
experiments under unknown interference, following Wang, Samii, Chang and
Aronow (2024).

## Usage

``` r
SpatialEffect(
  ras = NULL,
  Ydata = NULL,
  outcome = NULL,
  x_coord_Y = NULL,
  y_coord_Y = NULL,
  ras_Z = NULL,
  Zdata = NULL,
  x_coord_Z = NULL,
  y_coord_Z = NULL,
  treatment,
  covs = NULL,
  prob_treatment = NULL,
  dVec,
  dist.metric = "Euclidean",
  cType = "edge",
  numpts = NULL,
  evalpts = 1,
  only.unique = 0,
  smooth = 0,
  bw = NULL,
  bw_debias = NULL,
  bias_correction = TRUE,
  smooth.conley.se = 1,
  conf.band = 0,
  per.se = 1,
  blockvar = NULL,
  clustvar = NULL,
  conley.se = 1,
  kernel = "uni",
  cutoff = 0,
  alpha = 0.05,
  edf = FALSE,
  m = 2,
  lambda = 0.02,
  nPerms = 1000,
  n_threads = 1L,
  perm_engine = c("ri", "auto", "shuffle")
)
```

## Arguments

- ras:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of the outcome surface, or NULL to use kriging interpolation.

- Ydata:

  Data frame of outcome observations (used with kriging when outcome
  locations differ from intervention nodes).

- outcome:

  Character name of the outcome variable column.

- x_coord_Y, y_coord_Y:

  Character names of coordinate columns in Ydata.

- ras_Z:

  An [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object for
  polygon interventions, or NULL for point interventions. Can also be a
  SpatRaster or SpatVector which will be converted to sf.

- Zdata:

  Data frame of intervention node attributes.

- x_coord_Z, y_coord_Z:

  Character names of coordinate columns in Zdata (required for point
  interventions).

- treatment:

  Character name of the binary treatment column in Zdata.

- covs:

  Optional matrix of covariates (nz rows or nz\*length(dVec) rows).

- prob_treatment:

  Character name of propensity score column in Zdata for IPW
  (observational studies).

- dVec:

  Numeric vector of distance values at which to evaluate the AME.

- dist.metric:

  "Euclidean" (default) or "Geodesic".

- cType:

  Circle type: "edge" (default), "disk", or "donut".

- numpts:

  Number of sampling points on each circle (auto-selected if NULL).

- evalpts:

  Multiplier for auto-selected numpts (default 1).

- only.unique:

  Deduplicate raster cells per node (0 or 1).

- smooth:

  Whether to apply local polynomial smoothing (0 or 1).

- bw:

  Bandwidth for smoothing (auto-selected via CV if NULL).

- bw_debias:

  Bandwidth for bias correction (auto-selected if NULL).

- bias_correction:

  Use bias-corrected smoothing (default TRUE).

- smooth.conley.se:

  Compute Conley SE for smoothed estimates (0 or 1).

- conf.band:

  Compute uniform confidence bands for smoothed AME (0 or 1).

- per.se:

  Compute permutation-based CIs (0 or 1, default 1).

- blockvar:

  Block variable for stratified permutation.

- clustvar:

  Cluster variable for clustered permutation.

- conley.se:

  Compute Conley spatial HAC SEs (0 or 1, default 1).

- kernel:

  Kernel for HAC: "uni"/"uniform", "tri"/"triangular",
  "epa"/"epanechnikov".

- cutoff:

  Spatial bandwidth for Conley SE (0 = EHW).

- alpha:

  Significance level (default 0.05).

- edf:

  Use effective degrees of freedom adjustment (default FALSE).

- m:

  Polynomial order for kriging (default 2).

- lambda:

  Smoothing parameter for kriging (default 0.02).

- nPerms:

  Number of permutations (default 1000).

- n_threads:

  Number of threads for parallel C++ computation (default 1).

- perm_engine:

  Permutation generator for randomization inference: `"ri"` (default,
  uses [`ri::genperms()`](https://rdrr.io/pkg/ri/man/genperms.html)),
  `"auto"` (use `ri` if installed, otherwise fallback to shuffle), or
  `"shuffle"` (simple label reshuffling; ignores block/cluster design).

## Value

An S3 object of class "SpatialEffect" containing:

- AME_est:

  Data frame of distance and AME estimates

- Per.CI:

  Permutation confidence intervals (if per.se=1)

- Conley.SE:

  Conley standard errors (if conley.se=1)

- Conley.CI:

  Conley confidence intervals (if conley.se=1)

- AME_est_smoothed:

  Smoothed estimates (if smooth=1)

- smoothed.Conley.SE:

  Smoothed Conley SEs (if smooth=1 and smooth.conley.se=1)

- smoothed.Conley.CI:

  Smoothed Conley CIs

- smoothed.Conley.CB:

  Uniform confidence bands (if conf.band=1)

- Parameters:

  List of stored parameters for downstream use

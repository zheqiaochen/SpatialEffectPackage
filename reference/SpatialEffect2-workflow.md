# SpatialEffect2 Workflow Guide: Detailed Inputs, Parameters, and Practice

Practical operations manual for `SpatialEffect2`. This page documents
input templates, argument meanings, expected output formats, and
step-by-step analysis practice.

## Core functions

- [`SpatialEffect`](SpatialEffect.md): estimate distance-indexed AME
  curves.

- [`SpatialEffectTest`](SpatialEffectTest.md): randomization test on a
  distance range.

- S3 methods for class `"SpatialEffect"`:
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html),
  [`plot()`](https://rspatial.github.io/terra/reference/plot.html).

## Input templates (recommended formats)

**Template A: point interventions**

    # one row per intervention node
    Zdata <- data.frame(
      id = 1:nz,
      x = ...,   # numeric x coordinate
      y = ...,   # numeric y coordinate
      treat = ... # 0/1 treatment
    )

    # required call pieces
    x_coord_Z = "x"
    y_coord_Z = "y"
    treatment = "treat"

**Template B: polygon interventions**

    # ras_Z must be sf / SpatRaster / SpatVector
    ras_Z <- intervention_sf
    Zdata <- data.frame(treat = c(1,0,1,...))
    # node coordinates are derived internally from centroids

**Template C: kriging mode (no raster)**

    ras = NULL
    Ydata <- data.frame(x = ..., y = ..., outcome = ...)
    outcome = "outcome"
    x_coord_Y = "x"
    y_coord_Y = "y"

## Argument reference for `SpatialEffect`

**Always required**

- `Zdata`:

  Data frame of intervention nodes; one row per node.

- `treatment`:

  Character name of binary treatment column in `Zdata`.

- `dVec`:

  Numeric distance grid where AME is evaluated.

**Outcome and geometry inputs**

- `ras`:

  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  outcome surface. If `NULL`, kriging interpolation is used.

- `Ydata`:

  Outcome point data frame used in kriging mode.

- `outcome`:

  Outcome column name in `Ydata` (or in `Zdata` if `Ydata = NULL`).

- `x_coord_Y`, `y_coord_Y`:

  Coordinate column names in `Ydata`.

- `ras_Z`:

  Polygon interventions as `sf`, `SpatRaster`, or `SpatVector`; `NULL`
  means point intervention mode.

- `x_coord_Z`, `y_coord_Z`:

  Coordinate column names in `Zdata`; required in point mode.

**Distance definition and sampling**

- `cType`:

  Distance set type: `"edge"`, `"disk"`, or `"donut"`.

- `dist.metric`:

  `"Euclidean"` or `"Geodesic"`.

- `numpts`:

  Points sampled per circle boundary; auto-selected if `NULL`.

- `evalpts`:

  Multiplier for auto-selected `numpts`.

- `only.unique`:

  If `1`, deduplicate sampled raster cells per node.

**Design adjustment**

- `covs`:

  Optional covariate matrix. Accepts either `nz x p` or
  `(nz * length(dVec)) x p`.

- `prob_treatment`:

  Propensity-score column name for IPW adjustment.

- `blockvar`:

  Optional block variable for permutation design.

- `clustvar`:

  Optional cluster variable for permutation design.

**Inference controls**

- `per.se`:

  If `1`, compute permutation interval `Per.CI`.

- `conley.se`:

  If `1`, compute Conley HAC `Conley.SE` and `Conley.CI`.

- `kernel`:

  Conley kernel: `"uni"`/`"uniform"`, `"tri"`/`"triangular"`,
  `"epa"`/`"epanechnikov"`.

- `cutoff`:

  Spatial kernel bandwidth for Conley HAC.

- `alpha`:

  Significance level used for intervals.

- `edf`:

  Use effective-DOF adjustment for Conley SE.

- `nPerms`:

  Number of permutations.

**Smoothing controls**

- `smooth`:

  If `1`, run local-polynomial smoothing on AME curve.

- `bw`:

  Smoothing bandwidth; if `NULL`, selected by CV.

- `bw_debias`:

  Bias-correction bandwidth; if `NULL`, selected by CV when needed.

- `bias_correction`:

  If `TRUE`, compute debiased smoother.

- `smooth.conley.se`:

  If `1`, compute uncertainty for smoothed curve.

- `conf.band`:

  If `1`, compute uniform confidence band for smoothed curve.

**Kriging and compute controls**

- `m`, `lambda`:

  Kriging controls used when `ras = NULL`.

- `n_threads`:

  Thread count for C++ distance and HAC kernels.

## Argument reference for `SpatialEffectTest`

- `result.list`:

  A `"SpatialEffect"` object returned by `SpatialEffect`.

- `test.range`:

  Length-2 numeric vector `c(d_min, d_max)` for the tested distance
  interval.

- `smooth`:

  If `1`, test uses smoothed curve; otherwise unsmoothed AME.

- `alpha`:

  Significance level for permutation interval.

## Output format (example structure)

**`SpatialEffect(...)` returns class `"SpatialEffect"`.**

    $AMR_est               # data.frame with columns d, taud_est
    $Per.CI                # optional matrix [length(dVec), 2]
    $Conley.SE             # optional numeric vector [length(dVec)]
    $Conley.CI             # optional matrix [length(dVec), 2]
    $AMR_est_smoothed      # optional matrix [length(dVec), 2]
    $smoothed.Conley.SE    # optional numeric vector
    $smoothed.Conley.CI    # optional matrix [length(dVec), 2]
    $smoothed.Conley.CB    # optional matrix [length(dVec), 2]
    $Parameters            # list used by methods/tests
    $call                  # stored call

**`SpatialEffectTest(...)` returns:**

    $test.stat             # scalar statistic on selected range
    $test.CI               # length-2 permutation interval

## Practical workflow checklist

1.  Check CRS alignment and spatial overlap first.

2.  Start with conservative baseline: `smooth = 0`, `per.se = 1`,
    `conley.se = 1`.

3.  Use `dVec` step size similar to raster resolution for `"donut"` to
    reduce empty-ring `NA`.

4.  Add smoothing only after unsmoothed curve is interpretable.

5.  Test policy-relevant ranges with `SpatialEffectTest`.

6.  Run sensitivity checks over `cType`, `cutoff`, `kernel`, and `dVec`
    granularity.

## Common issues and fixes

**Issue: many `NA` values in AME**

- CRS mismatch between outcome and intervention objects.

- No or weak spatial overlap.

- Misalignment between `Zdata` rows and intervention geometries.

- `"donut"` rings too narrow relative to raster resolution.

**Issue: unstable smoothing / singular design**

- Increase `bw` and `bw_debias`.

- Reduce collinearity in `covs`.

- Use less dense `dVec`.

**Issue: runtime too slow**

- Increase `n_threads`.

- Decrease `nPerms`.

- Reduce sampling load via `numpts` and `evalpts`.

## Author

Package authors: Ye Wang, Cyrus Samii, Haoge Chang, and P. M. Aronow.

This workflow page is maintained as a practical user reference.

## Examples

``` r
if (FALSE) { # \dontrun{
library(SpatialEffect2)
library(terra)

set.seed(1)

# 1) outcome raster
ras <- rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 60)
values(ras) <- rnorm(ncell(ras))

# 2) intervention nodes
nz <- 100
Zdata <- data.frame(
  x = runif(nz, 2, 58),
  y = runif(nz, 2, 58),
  treat = rbinom(nz, 1, 0.5)
)

# 3) estimate AME
res <- SpatialEffect(
  ras = ras,
  Zdata = Zdata,
  x_coord_Z = "x",
  y_coord_Z = "y",
  treatment = "treat",
  dVec = seq(0, 10, by = 1),
  cType = "donut",
  dist.metric = "Euclidean",
  smooth = 1,
  per.se = 1,
  conley.se = 1,
  nPerms = 500,
  n_threads = 2L
)

summary(res, dVec.range = c(1, 5))
plot(res, ci.type = "both")
SpatialEffectTest(res, test.range = c(1, 5), smooth = 0)
} # }
```

# SpatialEffect2

`SpatialEffect2` estimates distance-indexed **Average Marginalized Effects (AME)** for spatial interventions under unknown interference, following Wang, Samii, Chang, and Aronow (2024).

It supports:
- Point or polygon interventions
- `edge` / `disk` / `donut` distance definitions
- Randomization-based uncertainty (`Per.CI`)
- Conley spatial HAC uncertainty (`Conley.CI`)
- Optional local-polynomial smoothing with bias correction

Tutorial is now written directly on the pkgdown home page:
- [pkgdown/index.md](pkgdown/index.md)

Build the pkgdown website from package root:
```r
pkgdown::build_site()
```

Chinese overview is available at:
- [README_zh.md](README_zh.md)

---

## Installation

From local source:

```r
install.packages("/path/to/SpatialEffect2", repos = NULL, type = "source")
```

Or from package root:

```bash
R CMD INSTALL SpatialEffect2
```

---

## API At A Glance

Core estimator:

```r
SpatialEffect(
  ras = NULL, Ydata = NULL, outcome = NULL,
  x_coord_Y = NULL, y_coord_Y = NULL,
  ras_Z = NULL, Zdata = NULL,
  x_coord_Z = NULL, y_coord_Z = NULL,
  treatment,
  covs = NULL, prob_treatment = NULL,
  dVec,
  dist.metric = "Euclidean",
  cType = "edge",
  numpts = NULL, evalpts = 1, only.unique = 0,
  smooth = 0, bw = NULL, bw_debias = NULL,
  bias_correction = TRUE, smooth.conley.se = 1, conf.band = 0,
  per.se = 1, blockvar = NULL, clustvar = NULL,
  conley.se = 1, kernel = "uni", cutoff = 0,
  alpha = 0.05, edf = FALSE,
  m = 2, lambda = 0.02,
  nPerms = 1000,
  n_threads = 1L
)
```

Randomization test for a policy-relevant distance range:

```r
SpatialEffectTest(result.list, test.range, smooth = 0, alpha = 0.05)
```

Display helpers:

```r
summary(result, dVec.range = c(1, 5))
plot(result, smooth = TRUE, ci.type = "both")
```

---

## Input Data Templates

### 1) Point interventions (`ras_Z = NULL`)

`Zdata` must contain at least:
- treatment assignment column (`0/1`)
- x and y coordinates of intervention nodes

Example:

```r
Zdata <- data.frame(
  id = 1:6,
  x = c(10, 15, 23, 30, 41, 44),
  y = c(12, 18, 20, 28, 33, 40),
  treat = c(1, 0, 1, 0, 1, 0)
)
```

Function call pieces:

```r
treatment = "treat",
x_coord_Z = "x",
y_coord_Z = "y"
```

### 2) Polygon interventions (`ras_Z` as `sf`)

If treatment units are polygons, pass polygon geometry in `ras_Z` and treatment in `Zdata` (same row order / row count). Node coordinates are computed from centroids internally.

```r
result <- SpatialEffect(
  ras = ras,
  ras_Z = intervention_sf,
  Zdata = data.frame(treat = c(1, 0, 1, 0)),
  treatment = "treat",
  dVec = seq(0, 20, by = 1),
  cType = "donut"
)
```

### 3) Kriging mode (`ras = NULL`)

When there is no full raster outcome surface, provide point outcomes:

```r
Ydata <- data.frame(
  x = runif(200, 0, 100),
  y = runif(200, 0, 100),
  outcome = rnorm(200)
)
```

Call pieces:

```r
ras = NULL,
Ydata = Ydata,
outcome = "outcome",
x_coord_Y = "x",
y_coord_Y = "y"
```

---

## Parameter Reference (`SpatialEffect`)

### A) Required arguments

- `Zdata` (`data.frame`): one row per intervention node.
- `treatment` (`character`): name of binary treatment column in `Zdata`.
- `dVec` (`numeric`): distance grid where AME is estimated.

Also required depending on geometry:
- Point intervention: `x_coord_Z`, `y_coord_Z` must be provided.
- Polygon intervention: provide `ras_Z` (`sf`, `SpatRaster`, or `SpatVector`).

### B) Outcome inputs

- `ras`: `terra::SpatRaster` outcome surface. If `NULL`, kriging is used.
- `Ydata`: outcome point table used when `ras = NULL`.
- `outcome`: outcome column name in `Ydata` (or in `Zdata` if `Ydata = NULL`).
- `x_coord_Y`, `y_coord_Y`: coordinate column names in `Ydata`.

### C) Distance and ring definition

- `cType`:
  - `"edge"`: boundary at distance `d`
  - `"disk"`: cumulative neighborhood `[0, d]`
  - `"donut"`: difference between adjacent cumulative rings
- `dist.metric`: `"Euclidean"` or `"Geodesic"`.
- `numpts`: number of points on each circle (auto if `NULL`).
- `evalpts`: multiplier for auto-selected `numpts`.
- `only.unique`: `1` deduplicates raster cells per node (`0` keeps repeats).

### D) Design adjustment and weighting

- `covs`: optional covariate matrix, either:
  - `nz x p` (will be repeated across `dVec`), or
  - `(nz * length(dVec)) x p` stacked by distance.
- `prob_treatment`: propensity-score column name in `Zdata` (IPW mode).
- `blockvar`, `clustvar`: optional design variables for permutations.

### E) Inference controls

- `per.se`: `1` computes permutation CIs (`Per.CI`).
- `conley.se`: `1` computes Conley SE/CI (`Conley.SE`, `Conley.CI`).
- `kernel`: Conley kernel, one of `"uni"`, `"tri"`, `"epa"` (or full names).
- `cutoff`: spatial kernel cutoff (`0` means EHW-like behavior in Conley routine).
- `alpha`: confidence level control (`0.05` gives 95% intervals).
- `edf`: effective degree-of-freedom adjustment for Conley SE.
- `nPerms`: number of permutations for randomization uncertainty.

### F) Smoothing controls

- `smooth`: `1` enables local-polynomial smoothing.
- `bw`: smoothing bandwidth (`NULL` triggers CV).
- `bw_debias`: bandwidth for bias correction (`NULL` triggers CV when needed).
- `bias_correction`: if `TRUE`, uses debiased smoother.
- `smooth.conley.se`: if `1`, compute uncertainty for smoothed curve.
- `conf.band`: if `1`, compute uniform confidence band for smoothed curve.

### G) Kriging and compute controls

- `m`, `lambda`: kriging parameters used when `ras = NULL`.
- `n_threads`: C++ thread count for distance and Conley routines.

---

## Output Format (What You Get Back)

`SpatialEffect(...)` returns an object of class `"SpatialEffect"`.

Main slots:

- `AME_est`: `data.frame` with columns:
  - `d`: distance value from `dVec`
  - `taud_est`: estimated AME at distance `d`
- `Per.CI` (optional): numeric matrix `length(dVec) x 2` (`lower`, `upper`).
- `Conley.SE` (optional): numeric vector `length(dVec)`.
- `Conley.CI` (optional): numeric matrix `length(dVec) x 2`.
- `AME_est_smoothed` (optional): numeric matrix with columns `(d, tau_smoothed)`.
- `smoothed.Conley.SE` / `smoothed.Conley.CI` / `smoothed.Conley.CB` (optional).
- `Parameters`: internal list used for downstream tests and methods.

Minimal output inspection:

```r
names(result)
# [1] "AME_est" "Per.CI" "Conley.SE" "Conley.CI" "AME_est_smoothed" ...

head(result$AME_est)
#     d   taud_est
# 1 0.0  0.012
# 2 1.0  0.019
# 3 2.0  0.024

dim(result$Conley.CI)
# [1] length(dVec) 2
```

`SpatialEffectTest(...)` returns:

- `test.stat`: scalar sum of AME over `test.range`.
- `test.CI`: length-2 permutation interval under the sharp null.

Example:

```r
test <- SpatialEffectTest(result, test.range = c(1, 5), smooth = 0)
str(test)
# List of 2
# $ test.stat: num ...
# $ test.CI  : Named num [1:2] ...
```

---

## Quick Start (Point Intervention)

```r
library(SpatialEffect2)
library(terra)

set.seed(1)

# Outcome raster
ras <- rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 60)
values(ras) <- rnorm(ncell(ras))

# Intervention nodes
nz <- 100
Zdata <- data.frame(
  x = runif(nz, 2, 58),
  y = runif(nz, 2, 58),
  treat = rbinom(nz, 1, 0.5)
)

result <- SpatialEffect(
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

summary(result, dVec.range = c(1, 5))
plot(result, ci.type = "both")

test <- SpatialEffectTest(result, test.range = c(1, 5), smooth = 0)
test
```

---

## Practical Reference: Recommended Workflow

1. Validate geometry first.
- Ensure CRS is consistent across outcome and intervention objects.
- Ensure intervention nodes/polygons spatially overlap outcome support.

2. Start with a conservative baseline.
- Use `cType = "donut"`, moderate `dVec` grid (for example, step near raster resolution), `per.se = 1`, `conley.se = 1`, `smooth = 0`.

3. Add smoothing only after baseline diagnostics.
- Turn on `smooth = 1` after confirming unsmoothed AME is stable.
- If CV bandwidths look unstable, set `bw` and `bw_debias` manually.

4. Define policy ranges and run tests.
- Use `SpatialEffectTest(result, test.range = c(a, b))` to evaluate practical distance ranges.

5. Run sensitivity checks.
- Compare `edge`/`disk`/`donut`.
- Re-run with alternative `cutoff`, `kernel`, and `dVec` granularity.
- Track whether conclusions change materially.

---

## Common Issues

### 1) AME is mostly `NA`

- Check CRS consistency.
- Check spatial overlap.
- Ensure `Zdata` rows match intervention geometry rows.
- For `donut`, increase `dVec` step or use finer outcome resolution.

### 2) Unstable smoothing or singular regressions

- Increase `bw` / `bw_debias`.
- Reduce collinear covariates.
- Use a less aggressive (less dense) distance grid.

### 3) Runtime is slow

- Increase `n_threads`.
- Reduce `nPerms`.
- Reduce circle sampling burden (`numpts`, `evalpts`).
- Fix bandwidths instead of CV if theory supports it.

---

## References

Wang, Y., Samii, C., Chang, H., and Aronow, P. M. (2024).  
Design-Based Inference for Spatial Experiments under Unknown Interference.

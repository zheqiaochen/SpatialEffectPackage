# SpatialEffect2 Tutorial


`SpatialEffect2` estimates distance-indexed **Average Marginalized Effects (AME)**
for spatial interventions under unknown interference.

This page is the main practical tutorial for the pkgdown site.

## 1. Setup

```r
library(SpatialEffect2)
library(terra)
```

## 2. Prepare Inputs

### Outcome surface (`ras`)

```r
set.seed(1)
ras <- rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 60)
values(ras) <- rnorm(ncell(ras))
```

### Point interventions (`Zdata`)

```r
nz <- 100
Zdata <- data.frame(
  x = runif(nz, 2, 58),
  y = runif(nz, 2, 58),
  treat = rbinom(nz, 1, 0.5)
)
```

Required columns for point mode:

- treatment column (`0/1`)
- x/y intervention coordinates

For polygon interventions, pass `ras_Z` as `sf`/`SpatVector`/`SpatRaster` and
keep one row in `Zdata` per polygon.

## 3. Choose Distance Design

- `cType = "edge"`: boundary at distance `d`
- `cType = "disk"`: cumulative region `[0, d]`
- `cType = "donut"`: adjacent-ring difference (often most interpretable)

```r
dVec <- seq(0, 10, by = 1)
```

Practical tip: under `donut`, choose `dVec` step size close to raster
resolution to reduce empty-ring `NA` issues.

## 4. Estimate AME Curve

```r
result <- SpatialEffect(
  ras = ras,
  Zdata = Zdata,
  x_coord_Z = "x",
  y_coord_Z = "y",
  treatment = "treat",
  dVec = dVec,
  cType = "donut",
  dist.metric = "Euclidean",
  smooth = 1,
  per.se = 1,
  conley.se = 1,
  kernel = "uni",
  cutoff = 0,
  alpha = 0.05,
  nPerms = 500,
  n_threads = 2L
)
```

## 5. Understand Outputs

`SpatialEffect(...)` returns class `"SpatialEffect"`.

Common slots:

- `AMR_est`: data frame with `d`, `taud_est`
- `Per.CI`: permutation CI matrix (`length(dVec) x 2`)
- `Conley.SE`, `Conley.CI`: spatial HAC uncertainty
- `AMR_est_smoothed`: smoothed AME curve
- `Parameters`: stored intermediate objects/settings

```r
summary(result, dVec.range = c(1, 5))
plot(result, ci.type = "both")
```

## 6. Test Policy-Relevant Distance Ranges

```r
test <- SpatialEffectTest(
  result.list = result,
  test.range = c(1, 5),
  smooth = 0,
  alpha = 0.05
)

test
```

## 7. Parameter Cheat Sheet

### Required

- `Zdata`
- `treatment`
- `dVec`

### Geometry and outcomes

- `ras` or (`Ydata`, `outcome`, `x_coord_Y`, `y_coord_Y`) for kriging mode
- `x_coord_Z`, `y_coord_Z` (point mode)
- `ras_Z` (polygon mode)

### Inference

- `per.se`, `nPerms`
- `conley.se`, `kernel`, `cutoff`, `edf`
- `alpha`

### Smoothing

- `smooth`, `bw`, `bw_debias`, `bias_correction`
- `smooth.conley.se`, `conf.band`

### Compute/performance

- `n_threads`, `numpts`, `evalpts`

## 8. Recommended Practice

1. Verify CRS consistency and spatial overlap before estimation.
2. Start with unsmoothed baseline (`smooth = 0`) and both uncertainty engines.
3. Add smoothing only after baseline AME is stable.
4. Test policy ranges with `SpatialEffectTest()`.
5. Run sensitivity checks over `dVec`, `cType`, `kernel`, and `cutoff`.

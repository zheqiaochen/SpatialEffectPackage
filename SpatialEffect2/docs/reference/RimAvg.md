# Compute circle averages for all intervention nodes at a given distance

Compute circle averages for all intervention nodes at a given distance

## Usage

``` r
RimAvg(
  ras,
  Yobs = NULL,
  ras_Z = NULL,
  nz,
  Zdata,
  x_coord_Z,
  y_coord_Z,
  dUp,
  numpts = NULL,
  gridRes = NULL,
  evalpts = 10,
  only.unique = 0,
  dtype = "raster",
  dist.metric = "Euclidean",
  cType = "edge",
  dist_unit = NULL,
  n_threads = 1L
)
```

## Arguments

- ras:

  SpatRaster or Krig object

- Yobs:

  Raster values vector (NULL for krig)

- ras_Z:

  sf object for polygon interventions (NULL for point)

- nz:

  Number of intervention nodes

- Zdata:

  Data frame of intervention data

- x_coord_Z:

  Name of x-coordinate column

- y_coord_Z:

  Name of y-coordinate column

- dUp:

  Distance value

- numpts:

  Number of sampling points

- gridRes:

  Raster resolution

- evalpts:

  Multiplier for numpts

- only.unique:

  Deduplicate grid cells?

- dtype:

  "raster" or "krig"

- dist.metric:

  "Euclidean" or "Geodesic"

- cType:

  "edge", "disk", or "donut"

- dist_unit:

  Step size for disk/donut

- n_threads:

  Number of threads for C++ computation

## Value

List(Ybard, Ybard_sum, Ybard_len)

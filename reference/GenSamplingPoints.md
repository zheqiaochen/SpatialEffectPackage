# Generate sampling points on a circle (edge), disk, or donut

Generate sampling points on a circle (edge), disk, or donut

## Usage

``` r
GenSamplingPoints(
  center,
  radius,
  numpts,
  dist.metric = "Euclidean",
  cType = "edge",
  dist_unit = NULL
)
```

## Arguments

- center:

  Numeric vector c(x, y)

- radius:

  Distance radius

- numpts:

  Number of points around the circle

- dist.metric:

  "Euclidean" or "Geodesic"

- cType:

  "edge", "disk", or "donut"

- dist_unit:

  Step size for disk/donut radii

## Value

Matrix of (x, y) sampling coordinates

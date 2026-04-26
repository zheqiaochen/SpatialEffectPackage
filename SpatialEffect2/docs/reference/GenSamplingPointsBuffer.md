# Generate sampling points by buffering a polygon (sf-based)

Generate sampling points by buffering a polygon (sf-based)

## Usage

``` r
GenSamplingPointsBuffer(geom, ras, radius, cType)
```

## Arguments

- geom:

  An sf geometry (single polygon)

- ras:

  A SpatRaster object

- radius:

  Buffer distance

- cType:

  "edge", "disk", or "donut"

## Value

Matrix of (x, y) sampling coordinates on valid raster cells

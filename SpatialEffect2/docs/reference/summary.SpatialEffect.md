# Summary method for SpatialEffect objects

Summarizes a `"SpatialEffect"` object as a compact table. The columns
are determined by the components available in the object: point
estimates are always shown, while Conley intervals, permutation
intervals, smoothed estimates, and smoothed intervals are shown only
when they were computed.

## Usage

``` r
# S3 method for class 'SpatialEffect'
summary(object, dVec.range = NULL, ...)
```

## Arguments

- object:

  A `"SpatialEffect"` object returned by
  [`SpatialEffect`](SpatialEffect.md).

- dVec.range:

  Optional numeric vector of length 2, `c(d_min, d_max)`. When supplied,
  only distances inside this closed interval are displayed. Use the same
  distance units as the original `dVec`.

- ...:

  Additional arguments passed by the generic. Currently ignored.

## Value

The formatted summary table is printed. The printed matrix is returned
invisibly by the underlying print call.

## See also

[`SpatialEffect`](SpatialEffect.md),
[`print.SpatialEffect`](print.SpatialEffect.md),
[`plot.SpatialEffect`](plot.SpatialEffect.md)

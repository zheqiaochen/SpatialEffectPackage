# Print method for SpatialEffect objects

Prints the distance-indexed AME estimates stored in a `"SpatialEffect"`
object. If uncertainty estimates were computed, the displayed table also
includes the available confidence-interval columns.

## Usage

``` r
# S3 method for class 'SpatialEffect'
print(x, dVec.range = NULL, ...)
```

## Arguments

- x:

  A `"SpatialEffect"` object returned by
  [`SpatialEffect`](SpatialEffect.md).

- dVec.range:

  Optional numeric vector of length 2, `c(d_min, d_max)`. When supplied,
  only distances inside this closed interval are displayed. Use the same
  distance units as the original `dVec`.

- ...:

  Additional arguments passed by the generic. Currently ignored.

## Value

The formatted results table is printed. The printed matrix is returned
invisibly by the underlying print call.

## See also

[`SpatialEffect`](SpatialEffect.md),
[`summary.SpatialEffect`](summary.SpatialEffect.md),
[`plot.SpatialEffect`](plot.SpatialEffect.md)

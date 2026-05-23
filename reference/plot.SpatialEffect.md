# Plot method for SpatialEffect objects

The default rendering (`style = "lines"`) is fully black-and-white-print
safe: confidence intervals are drawn as paired lines with distinct dash
patterns, so they remain legible in grayscale printing. Set
`style = "shade"` for translucent bands.

## Usage

``` r
# S3 method for class 'SpatialEffect'
plot(
  x,
  smooth = TRUE,
  ci.type = "conley",
  style = c("lines", "shade"),
  main = NULL,
  xlab = "Distance",
  ylab = "Average marginalized effect",
  ...
)
```

## Arguments

- x:

  A `"SpatialEffect"` object returned by
  [`SpatialEffect`](SpatialEffect.md).

- smooth:

  Logical. If `TRUE`, plot the smoothed AME curve when the object
  contains `AME_est_smoothed`. If no smoothed estimates are available,
  this argument is ignored.

- ci.type:

  Character string specifying which confidence intervals to draw:
  `"conley"`, `"permutation"`, or `"both"`. Requested intervals are
  drawn only if the corresponding components were computed by
  [`SpatialEffect`](SpatialEffect.md).

- style:

  Rendering style for confidence intervals: `"lines"` (default,
  B&W-print safe dashed lines) or `"shade"` (translucent polygons).

- main:

  Optional figure title. Defaults to no title (AER convention; titles
  typically appear as figure captions).

- xlab, ylab:

  Axis labels passed to the base plotting call.

- ...:

  Additional graphical arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `NULL`. The method is called for its plotting side
effect.

## See also

[`SpatialEffect`](SpatialEffect.md),
[`summary.SpatialEffect`](summary.SpatialEffect.md),
[`print.SpatialEffect`](print.SpatialEffect.md)

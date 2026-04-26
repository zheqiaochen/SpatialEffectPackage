# Plot method for SpatialEffect objects

Produces an AER-style figure. The default rendering (`style = "lines"`)
is fully black-and-white-print safe: confidence intervals are drawn as
paired lines with distinct dash patterns, so they remain legible in
grayscale printing. Set `style = "shade"` for translucent bands.

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

  A SpatialEffect object

- smooth:

  Plot smoothed estimates if available (default TRUE)

- ci.type:

  Which CI to plot: "conley", "permutation", or "both"

- style:

  Rendering style for confidence intervals: `"lines"` (default,
  B&W-print safe dashed lines) or `"shade"` (translucent polygons).

- main:

  Optional figure title. Defaults to no title (AER convention; titles
  typically appear as figure captions).

- xlab, ylab:

  Axis labels.

- ...:

  Additional arguments passed to plot()

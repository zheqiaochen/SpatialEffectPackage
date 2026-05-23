# SpatialEffect2 Analysis Workflow

A decision guide for planning a SpatialEffect2 analysis. The tutorial
shows how to run a complete example; this page helps decide what to
pass, which estimand to use, and which checks to run before interpreting
results.

## How to use this page

Use this page before or during an empirical analysis. It is not a second
tutorial. It is a checklist for choices that affect the estimand,
computation, and uncertainty.

## Decision 1: what is the outcome surface?

- Raster outcome:

  Use `ras` when the outcome is already a gridded surface, such as
  forest loss, pollution, temperature, or population density. This is
  the most direct workflow.

- Point outcome:

  Use `Ydata`, `outcome`, `x_coord_Y`, and `y_coord_Y` when outcomes are
  point measurements. Set `ras = NULL` so SpatialEffect2 interpolates an
  outcome surface.

**Check:** the outcome and intervention objects should use compatible
coordinate reference systems. For `dist.metric = "Euclidean"`, distances
are interpreted in the units of the projected CRS.

## Decision 2: what is the intervention geometry?

- Point interventions:

  Use `Zdata` with one row per node plus coordinate columns supplied
  through `x_coord_Z` and `y_coord_Z`.

- Polygon interventions:

  Pass polygons through `ras_Z`. Keep one row in `Zdata` per polygon, in
  the same order as `ras_Z`. Coordinates are taken from polygon
  centroids internally.

**Check:** `nrow(Zdata)` must match the number of intervention nodes.
For polygons, row order matters.

## Decision 3: which distance design answers the question?

- `cType = "edge"`:

  Use the outcome boundary at exactly distance `d`. This is closest to
  asking, "What happens right at this distance?"

- `cType = "disk"`:

  Use the cumulative region from `0` through `d`. This is useful for
  total neighborhood effects up to a radius.

- `cType = "donut"`:

  Use adjacent rings. This is often easiest to interpret when asking how
  effects change as we move outward.

**Check:** for `"donut"`, choose a `dVec` step that is not much smaller
than the raster resolution. Very thin rings can create empty-ring `NA`
values.

## Decision 4: experimental or observational study?

- Experimental study:

  If treatment probabilities are known by design and no weighting is
  needed, omit `prob_treatment`. Use `blockvar` or `clustvar` if
  randomization was blocked or clustered.

- Observational study:

  Estimate or supply a propensity score in `Zdata`, then pass its column
  name through `prob_treatment`. Inspect overlap before interpreting the
  weighted estimates.

**Check:** propensity scores near 0 or 1 can dominate weighted
estimates. Plot or summarize the propensity-score distribution.

## Decision 5: how should uncertainty be shown?

- `conley.se = 1`:

  Adds spatial HAC standard errors and `Conley.CI`. Use `cutoff` and
  `kernel` for sensitivity checks.

- `per.se = 1`:

  Adds permutation intervals `Per.CI`. Increase `nPerms` for final
  results.

- `smooth = 1`:

  Adds a smoothed AME curve. Use it after the unsmoothed curve is
  understandable.

- `smooth.conley.se = 1`:

  Adds uncertainty for the smoothed curve. This is more computationally
  demanding.

**Practical default:** start with `smooth = 0`, inspect the raw AME
curve, then add smoothing and sensitivity checks.

## Templates

**Point interventions**

    res <- SpatialEffect(
      ras = ras,
      Zdata = Zdata,
      x_coord_Z = "x",
      y_coord_Z = "y",
      treatment = "treat",
      dVec = dVec,
      cType = "donut"
    )

**Polygon interventions**

    res <- SpatialEffect(
      ras = ras,
      ras_Z = ras_Z,
      Zdata = Zdata,
      treatment = "treat",
      dVec = dVec,
      cType = "donut"
    )

**Observational study with propensity weights**

    res <- SpatialEffect(
      ras = ras,
      ras_Z = ras_Z,
      Zdata = Zdata,
      treatment = "treat",
      prob_treatment = "p_treat",
      dVec = dVec,
      cType = "donut"
    )

## After estimation

    summary(res, dVec.range = c(d_min, d_max))
    plot(res, ci.type = "conley")
    SpatialEffectTest(res, test.range = c(d_min, d_max), smooth = 0)

- [`summary()`](https://rdrr.io/pkg/terra/man/summary.html):

  Use for a compact numeric table at a relevant distance range.

- [`plot()`](https://rdrr.io/pkg/terra/man/plot.html):

  Use for the AME curve and available intervals.

- [`SpatialEffectTest()`](SpatialEffectTest.md):

  Use for a policy-relevant distance band, such as 0–5 km.

## Troubleshooting checklist

- Many `NA` estimates:

  Check CRS alignment, spatial overlap, row order between `ras_Z` and
  `Zdata`, and whether donut rings are too thin relative to raster
  resolution.

- Runtime is slow:

  Reduce `numpts`, reduce `nPerms`, start with `smooth = 0`, or increase
  `n_threads`.

- Confidence intervals are very wide:

  Check treatment balance, propensity-score overlap, distance bands with
  few outcome cells, and sensitivity to `cutoff`.

- Results change with `dVec`:

  Use a distance grid aligned with the raster resolution and report
  sensitivity to reasonable grids.

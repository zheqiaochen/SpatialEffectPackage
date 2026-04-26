# Sharp null randomization test for spatial treatment effects

Tests whether the AME is significantly different from zero over a
specified range of distances, using Fisher's sharp null hypothesis.

## Usage

``` r
SpatialEffectTest(result.list, test.range, smooth = 0, alpha = 0.05)
```

## Arguments

- result.list:

  A SpatialEffect object from [`SpatialEffect`](SpatialEffect.md).

- test.range:

  Numeric vector of length 2 specifying the distance range to test.

- smooth:

  Whether to use smoothed estimates (0 or 1).

- alpha:

  Significance level (default 0.05).

## Value

A list with:

- test.stat:

  The observed test statistic (sum of AME over test range)

- test.CI:

  Permutation confidence interval under the null

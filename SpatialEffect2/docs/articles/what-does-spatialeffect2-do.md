# What does SpatialEffect2 do?

`SpatialEffect2` answers one practical question:

> If treatment happens at one location, how do outcomes change nearby,
> farther away, and at each distance in between?

The answer is a curve. The horizontal axis is distance from the
intervention. The vertical axis is the estimated effect at that
distance.

## 1. The spatial spillover problem

Many interventions do not only affect the place where they are applied.

- A forest protection program may reduce logging inside a village
  boundary while moving pressure to nearby forests.
- A public health program in one village may affect behavior in
  neighboring villages.
- A policing program may move crime from one area to another area.

In these settings, an outcome at one place can depend on treatment at
many places. A simple treated-versus-control comparison can mix together
local effects, spillovers, displacement, and interactions among treated
locations.

`SpatialEffect2` separates this into a more direct question: how large
is the effect at distance `d` from an intervention location?

## 2. How SpatialEffect2 differs from common alternatives

### Parametric spatial models

One approach is to write down a spatial model, such as a spatial lag
model, autoregressive model, or model with a chosen spatial weight
matrix (Arbia 2006; Kelejian and Piras 2017).

That approach requires a strong choice about how effects move through
space. For example, the model may assume that effects decay smoothly
with distance. Real spillovers may be non-monotone: they may be strong
nearby, weaker in the middle, and stronger again farther away. Effects
from multiple treated locations may also add up or interact.

`SpatialEffect2` avoids choosing one spillover formula. It estimates an
effect curve directly from the observed outcome surface and the
treatment assignment.

### Exposure mapping approaches

Another approach is to create exposure categories, such as “inside a
treated polygon” or “within 5 km of any treated unit” (Aronow and Samii
2017).

This can work well when the exposure categories are known in advance. In
many applications, the right categories are unclear. A 5 km cutoff, a 10
km cutoff, and a “near any treated unit” rule can answer different
questions.

`SpatialEffect2` uses distance itself as the organizing variable.
Instead of choosing one exposure map, it estimates effects over a
sequence of distances.

## 3. Why this is useful

`SpatialEffect2` is useful because it gives a simple object to
interpret:

1.  **A curve over distance**: you can see whether effects are local,
    spread out, or displaced.
2.  **A clear target**: each point on the curve answers “what is the
    average effect at this distance?”
3.  **Less dependence on a spillover model**: you do not need to specify
    every pathway through which treatment spreads.
4.  **Built-in uncertainty tools**: the package can report permutation
    intervals and Conley spatial HAC intervals.

## 4. The core idea: average effects around each intervention

For each intervention location `i`, choose a distance `d`.

Then look at outcomes around that location:

- `edge`: average outcomes on the boundary at distance `d`;
- `disk`: average outcomes inside the region from `0` to `d`;
- `donut`: average outcomes in a ring, such as `[d - kappa, d]`.

This gives one “circle average” for node `i` at distance `d`.

Now compare two cases:

1.  node `i` is treated;
2.  node `i` is not treated.

Other nodes may also be treated. Instead of holding them fixed one by
one, `SpatialEffect2` averages over the treatment patterns that can
happen elsewhere. This is the key move. It makes the target effect
meaningful even when many locations can affect one another.

For node `i`, the effect at distance `d` is:

``` math
\tau_i(d) =
\text{average outcome around } i \text{ if } i \text{ is treated}
-
\text{average outcome around } i \text{ if } i \text{ is not treated}.
```

Then average this quantity over all intervention nodes:

``` math
\mathrm{AME}(d) = \frac{1}{N}\sum_{i=1}^N \tau_i(d).
```

This is the **Average Marginalized Effect** (AME). At `d = 0`, it
describes the effect at the intervention location. At larger `d`, it
describes how the effect changes as we move away from the intervention.

## 5. How estimation works in experiments

In a randomized experiment, treatment assignment gives a direct way to
compare treated and control intervention nodes.

For each distance `d`, `SpatialEffect2` does the following:

1.  Compute the circle average around every intervention node.
2.  Compare circle averages for treated and control nodes.
3.  Store the treatment-control difference as the AME estimate at
    distance `d`.
4.  Repeat this for every value in `dVec`.

In the package, this comparison is implemented as a regression of the
circle average on the treatment indicator. The treatment coefficient is
the estimated AME at that distance.

The package can also quantify uncertainty:

- `Per.CI`: randomization-based intervals from reassigning treatment
  labels;
- `Conley.SE` and `Conley.CI`: spatial HAC standard errors and intervals
  that allow nearby nodes to be statistically dependent;
- [`SpatialEffectTest()`](../reference/SpatialEffectTest.md): a
  randomization test for a distance range, such as `1` to `3` km.

## 6. How estimation works in observational studies

In observational studies, treatment is not randomized. Treated and
control locations may differ before treatment begins.

`SpatialEffect2` handles this through inverse-probability weights when
the user provides a propensity score with `prob_treatment`.

Let

``` math
e_i = P(Z_i = 1 \mid X_i)
```

be the probability that node `i` is treated, based on observed
covariates. The weight is

``` math
w_i = \frac{Z_i}{e_i} + \frac{1 - Z_i}{1 - e_i}.
```

These weights make treated and control nodes more comparable on the
observed covariates used to estimate `e_i`. After weighting, the package
estimates the same distance-by-distance AME curve.

The usual observational assumptions still matter:

- the important confounders are observed;
- treated and control nodes have enough overlap;
- the propensity scores are estimated well enough for the application.

## 7. A tiny worked example

### 7.1 AME as an average of node-level effects

``` r
set.seed(2026)
tau_i <- rnorm(8, mean = 0.03, sd = 0.01)
AME_hat <- mean(tau_i)

round(tau_i, 4)
#> [1] 0.0352 0.0192 0.0314 0.0292 0.0233 0.0048 0.0226 0.0198
round(AME_hat, 4)
#> [1] 0.0232
```

### 7.2 IPW weights used in observational settings

``` r
Z_i <- c(1, 0, 1, 0, 1)
e_i <- c(0.75, 0.35, 0.60, 0.40, 0.80)
w_i <- Z_i / e_i + (1 - Z_i) / (1 - e_i)

data.frame(
  unit = seq_along(Z_i),
  Z_i = Z_i,
  e_i = e_i,
  w_i = round(w_i, 4)
)
#>   unit Z_i  e_i    w_i
#> 1    1   1 0.75 1.3333
#> 2    2   0 0.35 1.5385
#> 3    3   1 0.60 1.6667
#> 4    4   0 0.40 1.6667
#> 5    5   1 0.80 1.2500
```

## 8. References

- Arbia, Giuseppe. 2006. *Spatial Econometrics: Statistical Foundations
  and Applications to Regional Convergence*.
- Aronow, Peter M., and Cyrus Samii. 2017. “Estimating Average Causal
  Effects Under General Interference, with Application to a Social
  Network Experiment.” *Annals of Applied Statistics* 11(4): 1912-1947.
- Hudgens, Michael G., and M. Elizabeth Halloran. 2008. “Toward Causal
  Inference with Interference.” *Journal of the American Statistical
  Association* 103(482): 832-842.
- Kelejian, Harry H., and Gianfranco Piras. 2017. *Spatial
  Econometrics*.
- Sävje, Fredrik, Peter M. Aronow, and Michael G. Hudgens. 2021.
  “Average Treatment Effects in the Presence of Unknown Interference.”
  *Annals of Statistics* 49(2): 673-701.
- Wang, Ye, Cyrus Samii, Haoge Chang, and P. M. Aronow. 2025.
  “Design-Based Inference for Spatial Experiments under Unknown
  Interference.” *The Annals of Appled Statistics*.

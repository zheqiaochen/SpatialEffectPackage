# Safe inverse of cross-product matrix

Uses a tiny ridge penalty fallback when `crossprod(X)` is singular or
nearly singular. If ridge escalation still fails, falls back to an
SVD-based pseudo-inverse.

## Usage

``` r
.safe_solve_crossprod(X, ridge_init = 1e-08, max_tries = 8L)
```

## Arguments

- X:

  Numeric matrix.

- ridge_init:

  Initial ridge value.

- max_tries:

  Maximum number of ridge escalation attempts.

## Value

A list with `inv` and `ridge`.

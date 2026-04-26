# Generate permutation matrix

Simple Bernoulli-style permutation matrix generator. Falls back to the
ri package if available and blockvar/clustvar are specified.

## Usage

``` r
.gen_perms(Zup, blockvar = NULL, clustvar = NULL, maxiter = 1000)
```

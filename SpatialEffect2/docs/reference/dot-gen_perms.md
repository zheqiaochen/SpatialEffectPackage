# Generate permutation matrix

Uses [`ri::genperms()`](https://rdrr.io/pkg/ri/man/genperms.html) by
default to preserve assignment mechanisms under complete, blocked, or
clustered randomization. A simple label reshuffle fallback is available
for users who explicitly request it.

## Usage

``` r
.gen_perms(
  Zup,
  blockvar = NULL,
  clustvar = NULL,
  maxiter = 1000,
  engine = c("ri", "auto", "shuffle")
)
```

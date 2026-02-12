# Update compromise or prospective/post-hoc power analysis without re-simulating

When a power or compromise analysis was performed in
[`Spower`](https://philchalmers.github.io/Spower/reference/Spower.md)
this function can be used to update the compromise or power criteria
without the need for re-simulating the experiment. For compromise
analyses a `beta_alpha` criteria must be supplied, while for
prospective/post-hoc power analyses the `sig.level` must be supplied.

## Usage

``` r
# S3 method for class 'Spower'
update(object, sig.level = 0.05, beta_alpha = NULL, predCI = 0.95, ...)
```

## Arguments

- object:

  object returned from
  [`Spower`](https://philchalmers.github.io/Spower/reference/Spower.md)
  where `power` was estimated or the `bete_alpha` criteria were supplied

- sig.level:

  Type I error rate (alpha)

- beta_alpha:

  Type II/Type I error ratio

- predCI:

  confidence interval precision (see
  [`Spower`](https://philchalmers.github.io/Spower/reference/Spower.md)
  for similar input)

- ...:

  arguments to be passed

## Value

object of class `Spower` with updated information

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

########
## Prospective power analysis update

# Estimate power using sig.level = .05 (default)
out <- p_t.test(n = 50, d = .5) |> Spower()

# update power estimate given sig.level=.01 and .20
update(out, sig.level=.01)
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n     d sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    50   0.5      0.01 NA   
#> 
#> Estimate of power: 0.459
#> 95% Confidence Interval: [0.446, 0.467]
update(out, sig.level=.20)
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n     d sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    50   0.5       0.2 NA   
#> 
#> Estimate of power: 0.880
#> 95% Confidence Interval: [0.876, 1.000]


########
## Compromise analysis update

# Solve beta/alpha ratio to specific error trade-off constant
out <- p_t.test(n = 50, d = .5) |> Spower(beta_alpha = 2)

# update beta_alpha criteria without re-simulating
update(out, beta_alpha=4)
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n     d sig.level power beta_alpha
#>   <dbl> <dbl>     <dbl> <lgl>      <dbl>
#> 1    50   0.5        NA NA             4
#> 
#> Estimate of Type I error rate (alpha/sig.level): 0.064
#> 95% Confidence Interval: [0.060, 0.069]
#> 
#> Estimate of power (1-beta): 0.742
#> 95% Confidence Interval: [0.734, 0.751]

# also works if compromise not initially run but prospective/post-hoc power was
out <- p_t.test(n = 50, d = .5) |> Spower()
update(out, beta_alpha=4)
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n     d sig.level power beta_alpha
#>   <dbl> <dbl>     <dbl> <lgl>      <dbl>
#> 1    50   0.5        NA NA             4
#> 
#> Estimate of Type I error rate (alpha/sig.level): 0.066
#> 95% Confidence Interval: [0.061, 0.071]
#> 
#> Estimate of power (1-beta): 0.737
#> 95% Confidence Interval: [0.728, 0.746]

# }
```

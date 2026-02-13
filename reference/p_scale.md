# p-value from Scale Test simulation

Simulates data given one or two parent distributions and returns a
p-value testing that the scale of the type distributions are the same.
Default implementation uses Gaussian distributions, however the
distribution function may be modified to reflect other populations of
interest. Uses [`ansari.test`](https://rdrr.io/r/stats/ansari.test.html)
or [`mood.test`](https://rdrr.io/r/stats/mood.test.html) for the
analysis.

## Usage

``` r
p_scale(
  n,
  scale,
  n2_n1 = 1,
  two.tailed = TRUE,
  exact = NULL,
  test = "Ansari",
  parent = function(n, ...) rnorm(n),
  ...,
  return_analysis = FALSE
)
```

## Arguments

- n:

  sample size per group

- scale:

  the scale to multiply the second group by (1 reflects equal scaling)

- n2_n1:

  sample size ratio

- two.tailed:

  logical; use two-tailed test?

- exact:

  a logical indicating whether an exact p-value should be computed

- test:

  type of method to use. Can be either `'Ansari'` or `'Mood'`

- parent:

  data generation function (default assumes Gaussian shape). Must be
  population mean centered

- ...:

  additional arguments to pass to simulation functions (if used)

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

## Value

a single p-value

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# n=30 per group,
#  Distributions Gaussian with sd=1 for first group and sd=2 for second
p_scale(30, scale=2)
#> [1] 3.110729e-05
p_scale(30, scale=2, test='Mood')
#> [1] 0.1173551

# compare chi-squared distributions
parent <- function(n, df, ...) rchisq(n, df=df) - df
p_scale(30, scale=2, parent=parent, df=3)
#> [1] 0.003431284

# \donttest{
  # empirical power of the experiments
  p_scale(30, scale=2) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:23
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n scale sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    30     2      0.05 NA   
#> 
#> Estimate of power: 0.774
#> 95% Confidence Interval: [0.765, 0.782]
  p_scale(30, scale=2, test='Mood') |> Spower()
#> 
#> Execution time (H:M:S): 00:00:02
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n scale test  sig.level power
#>   <dbl> <dbl> <chr>     <dbl> <lgl>
#> 1    30     2 Mood       0.05 NA   
#> 
#> Estimate of power: 0.867
#> 95% Confidence Interval: [0.860, 0.873]

  p_scale(30, scale=2, parent=parent, df=3) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:23
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n scale    df sig.level power
#>   <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1    30     2     3      0.05 NA   
#> 
#> Estimate of power: 0.907
#> 95% Confidence Interval: [0.901, 0.913]
  p_scale(30, scale=2, test='Mood', parent=parent, df=3) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:02
#> Design conditions: 
#> 
#> # A tibble: 1 × 6
#>       n scale test     df sig.level power
#>   <dbl> <dbl> <chr> <dbl>     <dbl> <lgl>
#> 1    30     2 Mood      3      0.05 NA   
#> 
#> Estimate of power: 0.949
#> 95% Confidence Interval: [0.945, 0.954]

# }
```

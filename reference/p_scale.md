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
#> [1] 0.04529126
p_scale(30, scale=2, test='Mood')
#> [1] 0.01834861

# compare chi-squared distributions
parent <- function(n, df, ...) rchisq(n, df=df) - df
p_scale(30, scale=2, parent=parent, df=3)
#> [1] 0.04212984

# \donttest{
  # empirical power of the experiments
  p_scale(30, scale=2) |> Spower()
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 4
#>       n scale sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    30     2      0.05 NA   
#> 
#> Estimate of power: 0.776
#> 95% Confidence Interval: [0.768, 0.785]
#> Execution time (H:M:S): 00:00:21
  p_scale(30, scale=2, test='Mood') |> Spower()
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 5
#>       n scale test  sig.level power
#>   <dbl> <dbl> <chr>     <dbl> <lgl>
#> 1    30     2 Mood       0.05 NA   
#> 
#> Estimate of power: 0.866
#> 95% Confidence Interval: [0.859, 0.872]
#> Execution time (H:M:S): 00:00:02

  p_scale(30, scale=2, parent=parent, df=3) |> Spower()
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 5
#>       n scale    df sig.level power
#>   <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1    30     2     3      0.05 NA   
#> 
#> Estimate of power: 0.904
#> 95% Confidence Interval: [0.899, 0.910]
#> Execution time (H:M:S): 00:00:22
  p_scale(30, scale=2, test='Mood', parent=parent, df=3) |> Spower()
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 6
#>       n scale test     df sig.level power
#>   <dbl> <dbl> <chr> <dbl>     <dbl> <lgl>
#> 1    30     2 Mood      3      0.05 NA   
#> 
#> Estimate of power: 0.945
#> 95% Confidence Interval: [0.941, 0.950]
#> Execution time (H:M:S): 00:00:02

# }
```

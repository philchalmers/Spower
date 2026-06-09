# p-value from Shapiro-Wilk Normality Test simulation

Generates univariate distributional data and returns a p-value to assess
the null that the population follows a Gaussian distribution shape. Uses
[`shapiro.test`](https://rdrr.io/r/stats/shapiro.test.html).

## Usage

``` r
p_shapiro.test(dist, return_analysis = FALSE)
```

## Arguments

- dist:

  expression used to generate the required sample data

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

## Value

a single p-value

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

# 50 observations drawn from normal distribution (null is true)
p_shapiro.test(rnorm(50))
#> [1] 0.647198

# return analysis object
p_shapiro.test(rnorm(50), TRUE)
#> 
#>  Shapiro-Wilk normality test
#> 
#> data:  dist
#> W = 0.96395, p-value = 0.1301
#> 

# 50 observations from slightly skewed chi-squared distribution (power)
p_shapiro.test(rchisq(50, df=100))
#> [1] 0.1067649

# \donttest{
    # empirical Type I error rate estimate
    p_shapiro.test(rnorm(50)) |> Spower()
#> Warning: number of items to replace is not a multiple of replacement length
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 3
#>       dist sig.level power
#>      <dbl>     <dbl> <lgl>
#> 1 -0.27804      0.05 NA   
#> 
#> Estimate of power: 0.052
#> 95% Confidence Interval: [0.048, 0.057]
#> Execution time (H:M:S): 00:00:01

    # power
    p_shapiro.test(rchisq(50, df=100)) |> Spower()
#> Warning: number of items to replace is not a multiple of replacement length
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 3
#>     dist sig.level power
#>    <dbl>     <dbl> <lgl>
#> 1 96.174      0.05 NA   
#> 
#> Estimate of power: 0.108
#> 95% Confidence Interval: [0.102, 0.114]
#> Execution time (H:M:S): 00:00:01
# }
```

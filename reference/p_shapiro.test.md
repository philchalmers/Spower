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
#> [1] 0.624697

# return analysis object
p_shapiro.test(rnorm(50), TRUE)
#> 
#>  Shapiro-Wilk normality test
#> 
#> data:  dist
#> W = 0.95519, p-value = 0.05597
#> 

# 50 observations from slightly skewed chi-squared distribution (power)
p_shapiro.test(rchisq(50, df=100))
#> [1] 0.9646621

# \donttest{
    # empirical Type I error rate estimate
    p_shapiro.test(rnorm(50)) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:02
#> Design conditions: 
#> 
#> # A tibble: 1 × 2
#>   sig.level power
#>       <dbl> <lgl>
#> 1      0.05 NA   
#> 
#> Estimate of power: 0.051
#> 95% Confidence Interval: [0.047, 0.056]

    # power
    p_shapiro.test(rchisq(50, df=100)) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:02
#> Design conditions: 
#> 
#> # A tibble: 1 × 2
#>   sig.level power
#>       <dbl> <lgl>
#> 1      0.05 NA   
#> 
#> Estimate of power: 0.102
#> 95% Confidence Interval: [0.096, 0.108]
# }
```

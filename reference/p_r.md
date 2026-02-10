# p-value from correlation simulation

Generates correlated X-Y data and returns a p-value to assess the null
of no correlation in the population. The X-Y data are generated assuming
a bivariate normal distribution.

## Usage

``` r
p_r(
  n,
  r,
  rho = 0,
  method = "pearson",
  two.tailed = TRUE,
  gen_fun = gen_r,
  return_analysis = FALSE,
  ...
)

gen_r(n, r, ...)
```

## Arguments

- n:

  sample size

- r:

  correlation

- rho:

  population coefficient to test against. Uses the Fisher's
  z-transformation approximation when non-zero

- method:

  method to use to compute the correlation (see
  [`cor.test`](https://rdrr.io/r/stats/cor.test.html)). Only used when
  `rho = 0`

- two.tailed:

  logical; should a two-tailed or one-tailed test be used?

- gen_fun:

  function used to generate the required dependent bivariate data.
  Object returned must be a `matrix` with two columns and `n` rows.
  Default uses `gen_r` to generate conditionally dependent data from a
  bivariate normal distribution. User defined version of this function
  must include the argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization? Note that if `rho != 0` the `p.value` and related
  element will be replaced with internally computed approximation
  versions

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

## Value

a single p-value

## See also

`gen_r`

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# 50 observations, .5 correlation
p_r(50, r=.5)
#> [1] 1.118634e-05
p_r(50, r=.5, method = 'spearman')
#> [1] 8.326385e-05

# test against constant other than rho = .6
p_r(50, .5, rho=.60)
#> [1] 0.2973482

# return analysis model
p_r(50, .5, return_analysis=TRUE)
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  x and y
#> t = 3.2934, df = 48, p-value = 0.001863
#> alternative hypothesis: true correlation is not equal to 0
#> 95 percent confidence interval:
#>  0.1714634 0.6321291
#> sample estimates:
#>       cor 
#> 0.4293211 
#> 
p_r(50, .5, rho=.60, return_analysis=TRUE)
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  x and y
#> t = -0.38724, df = Inf, p-value = 0.6986
#> alternative hypothesis: true correlation is not equal to 0.6
#> 95 percent confidence interval:
#>  0.3370604 0.7271028
#> sample estimates:
#>       cor 
#> 0.5626228 
#> 

# \donttest{
    # compare simulated results to pwr package

    pwr::pwr.r.test(r=0.3, n=50)
#> 
#>      approximate correlation power calculation (arctangh transformation) 
#> 
#>               n = 50
#>               r = 0.3
#>       sig.level = 0.05
#>           power = 0.5715558
#>     alternative = two.sided
#> 
    p_r(n=50, r=0.3) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:08
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n     r sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    50   0.3      0.05 NA   
#> 
#> Estimate of power: 0.575
#> 95% Confidence Interval: [0.566, 0.585]

    pwr::pwr.r.test(r=0.3, power=0.80)
#> 
#>      approximate correlation power calculation (arctangh transformation) 
#> 
#>               n = 84.07364
#>               r = 0.3
#>       sig.level = 0.05
#>           power = 0.8
#>     alternative = two.sided
#> 
    p_r(n=NA, r=0.3) |> Spower(power=.80, interval=c(10, 200))
#> 
#> Execution time (H:M:S): 00:00:32
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n     r sig.level power
#>   <dbl> <dbl>     <dbl> <dbl>
#> 1    NA   0.3      0.05   0.8
#> 
#> Estimate of n: 84.0
#> 95% Predicted Confidence Interval: [83.1, 84.8]

    pwr::pwr.r.test(r=0.1, power=0.80)
#> 
#>      approximate correlation power calculation (arctangh transformation) 
#> 
#>               n = 781.7516
#>               r = 0.1
#>       sig.level = 0.05
#>           power = 0.8
#>     alternative = two.sided
#> 
    p_r(n=NA, r=0.1) |> Spower(power=.80, interval=c(200, 1000))
#> 
#> Execution time (H:M:S): 00:01:16
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n     r sig.level power
#>   <dbl> <dbl>     <dbl> <dbl>
#> 1    NA   0.1      0.05   0.8
#> 
#> Estimate of n: 788.1
#> 95% Predicted Confidence Interval: [772.8, 802.1]

# }
```

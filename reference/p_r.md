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

gen_r(n, r, means = c(0, 0), ...)
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

- means:

  mean vector for data generation, of length 2. Defaults to `c(0,0)`

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
#> [1] 1.863372e-05
p_r(50, r=.5, method = 'spearman')
#> [1] 0.01228988

# test against constant other than rho = .6
p_r(50, .5, rho=.60)
#> [1] 0.1723153

# return analysis model
p_r(50, .5, return_analysis=TRUE)
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  x and y
#> t = 3.8027, df = 48, p-value = 0.0004044
#> alternative hypothesis: true correlation is not equal to 0
#> 95 percent confidence interval:
#>  0.2341751 0.6698012
#> sample estimates:
#>       cor 
#> 0.4811598 
#> 
p_r(50, .5, rho=.60, return_analysis=TRUE)
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  x and y
#> t = -1.0488, df = Inf, p-value = 0.2942
#> alternative hypothesis: true correlation is not equal to 0.6
#> 95 percent confidence interval:
#>  0.2489252 0.6783477
#> sample estimates:
#>       cor 
#> 0.4931067 
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
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 4
#>       n     r sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    50   0.3      0.05 NA   
#> 
#> Estimate of power: 0.571
#> 95% Confidence Interval: [0.561, 0.581]
#> Execution time (H:M:S): 00:00:07

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
    p_r(n=interval(10, 200), r=0.3) |> Spower(power=.80)
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 4
#>       n     r sig.level power
#>   <dbl> <dbl>     <dbl> <dbl>
#> 1    NA   0.3      0.05   0.8
#> 
#> Estimate of n: 83.9
#> 95% Confidence Interval: [82.6, 85.2]
#> Execution time (H:M:S): 00:00:55

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
    p_r(n=interval(200, 1000), r=0.1) |> Spower(power=.80)
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 4
#>       n     r sig.level power
#>   <dbl> <dbl>     <dbl> <dbl>
#> 1    NA   0.1      0.05   0.8
#> 
#> Estimate of n: 781.8
#> 95% Confidence Interval: [775.1, 788.1]
#> Execution time (H:M:S): 00:00:34

# }
```

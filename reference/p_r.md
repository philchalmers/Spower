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
#> [1] 6.934215e-05
p_r(50, r=.5, method = 'spearman')
#> [1] 3.08572e-05

# test against constant other than rho = .6
p_r(50, .5, rho=.60)
#> [1] 0.0006483459

# return analysis model
p_r(50, .5, return_analysis=TRUE)
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  x and y
#> t = 3.3509, df = 48, p-value = 0.001576
#> alternative hypothesis: true correlation is not equal to 0
#> 95 percent confidence interval:
#>  0.1787245 0.6366053
#> sample estimates:
#>       cor 
#> 0.4354115 
#> 
p_r(50, .5, rho=.60, return_analysis=TRUE)
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  x and y
#> t = 0.83644, df = Inf, p-value = 0.4029
#> alternative hypothesis: true correlation is not equal to 0.6
#> 95 percent confidence interval:
#>  0.4848187 0.8008740
#> sample estimates:
#>      cor 
#> 0.672424 
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
#> Estimate of power: 0.580
#> 95% Confidence Interval: [0.570, 0.590]
#> Execution time (H:M:S): 00:00:08

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
#> Estimate of n: 83.7
#> 95% Confidence Interval: [83.1, 84.4]
#> Execution time (H:M:S): 00:00:32

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
#> Estimate of n: 774.7
#> 95% Confidence Interval: [768.1, 782.0]
#> Execution time (H:M:S): 00:00:36

# }
```

# p-value from variance test simulation

Generates one or or more sets of continuous data group-level data to
perform a variance test, and return a p-value. When two-samples are
investigated the [`var.test`](https://rdrr.io/r/stats/var.test.html)
function will be used, otherwise functions from the `EnvStats` package
will be used.

## Usage

``` r
p_var.test(
  n,
  vars,
  n.ratios = rep(1, length(vars)),
  sigma2 = 1,
  two.tailed = TRUE,
  test = "Levene",
  correct = TRUE,
  gen_fun = gen_var.test,
  return_analysis = FALSE,
  ...
)

gen_var.test(n, vars, n.ratios = rep(1, length(vars)), ...)
```

## Arguments

- n:

  sample size per group, assumed equal across groups

- vars:

  a vector of variances to use for each group; length of 1 for
  one-sample tests

- n.ratios:

  allocation ratios reflecting the sample size ratios. Default of 1 sets
  the groups to be the same size (n \* n.ratio)

- sigma2:

  population variance to test against in one-sample test

- two.tailed:

  logical; should a two-tailed or one-tailed test be used?

- test:

  type of test to use in multi-sample applications. Can be either
  `'Levene'` (default), `'Bartlett'`, or `'Fligner'`

- correct:

  logical; use correction when `test = 'Bartlett'`?

- gen_fun:

  function used to generate the required discrete data. Object returned
  must be a `matrix` with k rows and k columns of counts. Default uses
  `gen_var.test`. User defined version of this function must include the
  argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

## Value

a single p-value

## See also

`gen_var.test`

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# one sample
p_var.test(100, vars=10, sigma2=9)
#> [1] 0.8698271

# return analysis object
p_var.test(100, vars=10, sigma2=9, return_analysis = TRUE)
#> $statistic
#> Chi-Squared 
#>    98.47016 
#> 
#> $parameters
#> df 
#> 99 
#> 
#> $p.value
#> [1] 0.9922601
#> 
#> $estimate
#> variance 
#> 8.951833 
#> 
#> $null.value
#> variance 
#>        9 
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $method
#> [1] "Chi-Squared Test on Variance"
#> 
#> $data.name
#> [1] "dat$DV"
#> 
#> $conf.int
#>       LCL       UCL 
#>  6.900932 12.080404 
#> attr(,"conf.level")
#> [1] 0.95
#> 
#> attr(,"class")
#> [1] "htestEnvStats"

# three sample
p_var.test(100, vars=c(10, 9, 11))
#> [1] 0.1659917
p_var.test(100, vars=c(10, 9, 11), test = 'Fligner')
#> [1] 0.3943148
p_var.test(100, vars=c(10, 9, 11), test = 'Bartlett')
#> [1] 0.6316487

# \donttest{
  # power to detect three-group variance differences
  p_var.test(n=100, vars=c(10,9,11)) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:34
#> Design conditions: 
#> 
#> # A tibble: 1 × 3
#>       n sig.level power
#>   <dbl>     <dbl> <lgl>
#> 1   100      0.05 NA   
#> 
#> Estimate of power: 0.119
#> 95% Confidence Interval: [0.113, 0.125]

  # sample size per group to achieve 80% power
  p_var.test(n=NA, vars=c(10,9,11)) |>
         Spower(power=.80, interval=c(100, 1000))
#> 
#> Execution time (H:M:S): 00:04:41
#> Design conditions: 
#> 
#> # A tibble: 1 × 3
#>       n sig.level power
#>   <dbl>     <dbl> <dbl>
#> 1    NA      0.05   0.8
#> 
#> Estimate of n: 999.0
#> 95% Predicted Confidence Interval: [NA, NA]
# }
```

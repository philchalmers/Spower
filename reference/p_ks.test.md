# p-value from Kolmogorov-Smirnov one- or two-sample simulation

Generates one or two sets of continuous data group-level data and
returns a p-value under the null that the groups were drawn from the
same distribution (two sample) or from a theoretically known
distribution (one sample).

## Usage

``` r
p_ks.test(
  n,
  p1,
  p2,
  n2_n1 = 1,
  two.tailed = TRUE,
  parent = NULL,
  return_analysis = FALSE,
  ...
)
```

## Arguments

- n:

  sample size per group, assumed equal across groups

- p1:

  a function indicating how the data were generated for group 1

- p2:

  (optional) a function indicating how the data were generated for
  group 2. If omitted a one-sample test will be evaluated provided that
  `parent` is also specified

- n2_n1:

  sample size ratio. Default uses equal sample sizes

- two.tailed:

  logical; should a two-tailed or one-tailed test be used?

- parent:

  the cumulative distribution function to use (e.g.,
  [`pnorm`](https://rdrr.io/r/stats/Normal.html)). Specifying this input
  will construct a one-sample test setup

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to the `parent` distribution
  function from [`ks.test`](https://rdrr.io/r/stats/ks.test.html), as
  well as any other relevant parameter to `ks.test` (e.g.,
  `exact = TRUE`)

## Value

a single p-value

## See also

[`gen_t.test`](https://philchalmers.github.io/Spower/reference/p_t.test.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# two-sample test from two Gaussian distributions with different locations
p1 <- function(n) rnorm(n)
p2 <- function(n) rnorm(n, mean=-.5)
p_ks.test(n=100, p1, p2)
#> [1] 0.02431031

# return analysis model
p_ks.test(n=100, p1, p2, return_analysis=TRUE)
#> 
#>  Asymptotic two-sample Kolmogorov-Smirnov test
#> 
#> data:  dat1 and dat2
#> D = 0.17, p-value = 0.1111
#> alternative hypothesis: two-sided
#> 

# one-sample data from chi-squared distribution tested
#   against a standard normal distribution
pc <- function(n, df=15) (rchisq(n, df=df) - df) / sqrt(2*df)
p_ks.test(n=100, p1=pc, parent=pnorm, mean=0, sd=1)
#> [1] 0.05754102

# \donttest{
  # empirical power estimates
  p_ks.test(n=100, p1, p2) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 3
#>       n sig.level power
#>   <dbl>     <dbl> <lgl>
#> 1   100      0.05 NA   
#> 
#> Estimate of power: 0.828
#> 95% Confidence Interval: [0.821, 0.835]
  p_ks.test(n=100, p1=pc, parent=pnorm, mean=0, sd=1) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n  mean    sd sig.level power
#>   <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1   100     0     1      0.05 NA   
#> 
#> Estimate of power: 0.144
#> 95% Confidence Interval: [0.137, 0.151]

# }
```

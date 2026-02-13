# p-value from chi-squared test simulation

Generates multinomial data suitable for analysis with
[`chisq.test`](https://rdrr.io/r/stats/chisq.test.html).

## Usage

``` r
p_chisq.test(
  n,
  w,
  df,
  correct = TRUE,
  P0 = NULL,
  P = NULL,
  gen_fun = gen_chisq.test,
  return_analysis = FALSE,
  ...
)

gen_chisq.test(n, P, ...)
```

## Arguments

- n:

  sample size per group

- w:

  Cohen's w effect size

- df:

  degrees of freedom

- correct:

  logical; apply continuity correction?

- P0:

  specific null pattern, specified as a numeric vector or matrix

- P:

  specific power configuration, specified as a numeric vector or matrix

- gen_fun:

  function used to generate the required discrete data. Object returned
  must be a `matrix` with k rows and k columns of counts. Default uses
  `gen_chisq.test`. User defined version of this function must include
  the argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

## Value

a single p-value

## See also

`gen_chisq.test`

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# effect size w + df
p_chisq.test(100, w=.2, df=3)
#> [1] 0.6776621

# return analysis model
p_chisq.test(100, w=.2, df=3, return_analysis=TRUE)
#> 
#>  Chi-squared test for given probabilities
#> 
#> data:  tab
#> X-squared = 2.48, df = 3, p-value = 0.4789
#> 

# vector of explicit probabilities (goodness of fit test)
p_chisq.test(100, P0 = c(.25, .25, .25, .25),
                   P = c(.6, .2, .1, .1))
#> [1] 6.683395e-20

# matrix of explicit probabilities (two-dimensional test of independence)
p_chisq.test(100, P0 = matrix(c(.25, .25, .25, .25), 2, 2),
                   P = matrix(c(.6, .2, .1, .1),2,2))
#> Warning: Chi-squared approximation may be incorrect
#> [1] 0.1724933

# \donttest{
    # compare simulated results to pwr package

    P0 <- c(1/3, 1/3, 1/3)
    P <- c(.5, .25, .25)
    w <- pwr::ES.w1(P0, P)
    df <- 3-1
    pwr::pwr.chisq.test(w=w, df=df, N=100, sig.level=0.05)
#> 
#>      Chi squared power calculation 
#> 
#>               w = 0.3535534
#>               N = 100
#>              df = 2
#>       sig.level = 0.05
#>           power = 0.8962428
#> 
#> NOTE: N is the number of observations
#> 

    # slightly less power when evaluated empirically
    p_chisq.test(n=100, w=w, df=df) |> Spower(replications=100000)
#> 
#> Execution time (H:M:S): 00:00:22
#> Design conditions: 
#> 
#> # A tibble: 1 × 3
#>       n sig.level power
#>   <dbl>     <dbl> <lgl>
#> 1   100      0.05 NA   
#> 
#> Estimate of power: 0.886
#> 95% Confidence Interval: [0.884, 0.888]
    p_chisq.test(n=100, P0=P0, P=P) |> Spower(replications=100000)
#> 
#> Execution time (H:M:S): 00:00:17
#> Design conditions: 
#> 
#> # A tibble: 1 × 3
#>       n sig.level power
#>   <dbl>     <dbl> <lgl>
#> 1   100      0.05 NA   
#> 
#> Estimate of power: 0.886
#> 95% Confidence Interval: [0.884, 0.888]

    # slightly differ (latter more conservative due to finite sampling behaviour)
    pwr::pwr.chisq.test(w=w, df=df, power=.8, sig.level=0.05)
#> 
#>      Chi squared power calculation 
#> 
#>               w = 0.3535534
#>               N = 77.07751
#>              df = 2
#>       sig.level = 0.05
#>           power = 0.8
#> 
#> NOTE: N is the number of observations
#> 
    p_chisq.test(n=interval(50, 200), w=w, df=df) |> Spower(power=.80)
#> 
#> Execution time (H:M:S): 00:00:21
#> Design conditions: 
#> 
#> # A tibble: 1 × 3
#>       n sig.level power
#>   <dbl>     <dbl> <dbl>
#> 1    NA      0.05   0.8
#> 
#> Estimate of n: 79.9
#> 95% Predicted Confidence Interval: [79.4, 80.5]
    p_chisq.test(n=interval(50, 200), w=w, df=df, correct=FALSE) |>
      Spower(power=.80)
#> 
#> Execution time (H:M:S): 00:00:16
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n correct sig.level power
#>   <dbl> <lgl>       <dbl> <dbl>
#> 1    NA FALSE        0.05   0.8
#> 
#> Estimate of n: 80.1
#> 95% Predicted Confidence Interval: [79.3, 80.9]

    # Spower slightly more conservative even with larger N
    pwr::pwr.chisq.test(w=.1, df=df, power=.95, sig.level=0.05)
#> 
#>      Chi squared power calculation 
#> 
#>               w = 0.1
#>               N = 1544.324
#>              df = 2
#>       sig.level = 0.05
#>           power = 0.95
#> 
#> NOTE: N is the number of observations
#> 
    p_chisq.test(n=interval(1000, 2000), w=.1, df=df) |> Spower(power=.95)
#> 
#> Execution time (H:M:S): 00:00:08
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n     w sig.level power
#>   <dbl> <dbl>     <dbl> <dbl>
#> 1    NA   0.1      0.05  0.95
#> 
#> Estimate of n: 1568.5
#> 95% Predicted Confidence Interval: [1543.8, 1593.3]
    p_chisq.test(n=interval(1000, 2000), w=.1, df=df, correct=FALSE) |>
           Spower(power=.95)
#> 
#> Execution time (H:M:S): 00:00:07
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n     w correct sig.level power
#>   <dbl> <dbl> <lgl>       <dbl> <dbl>
#> 1    NA   0.1 FALSE        0.05  0.95
#> 
#> Estimate of n: 1574.5
#> 95% Predicted Confidence Interval: [1570.5, 1579.7]

# }
```

# p-value from three-variable mediation analysis simulation

Simple 3-variable mediation analysis simulation to test the hypothesis
that X -\> Y is mediated by the relationship X -\> M -\> Y. Currently, M
and Y are assumed to be continuous variables with Gaussian errors, while
X may be continuous or dichotomous.

## Usage

``` r
p_mediation(
  n,
  a,
  b,
  cprime,
  dichotomous.X = FALSE,
  two.tailed = TRUE,
  method = "wald",
  sd.X = 1,
  sd.Y = 1,
  sd.M = 1,
  gen_fun = gen_mediation,
  return_analysis = FALSE,
  ...
)

gen_mediation(
  n,
  a,
  b,
  cprime,
  dichotomous.X = FALSE,
  sd.X = 1,
  sd.Y = 1,
  sd.M = 1,
  ...
)
```

## Arguments

- n:

  total sample size unless `dichotomous.X = TRUE`, in which the value
  represents the size per group

- a:

  regression coefficient for the path X -\> M

- b:

  regression coefficient for the path M -\> Y

- cprime:

  partial regression coefficient for the path X -\> Y

- dichotomous.X:

  logical; should the X variable be generated as though it were
  dichotomous? If TRUE then `n` represents the sample size per group

- method:

  type of inferential method to use. Default uses the Wald (a.k.a.,
  Sobel) test

- sd.X:

  standard deviation for X

- sd.Y:

  standard deviation for Y

- sd.M:

  standard deviation for M

- gen_fun:

  function used to generate the required two-sample data. Object
  returned must be a `data.frame` with the columns `"DV"` and `"group"`.
  Default uses `gen_mediation` to generate conditionally Gaussian
  distributed samples. User defined version of this function must
  include the argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

## Value

a single p-value

## See also

`gen_mediation`

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

# joint test H0: a*b = 0
p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39)
#> [1] 1.012941e-05
p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, dichotomous.X=TRUE)
#> [1] 1.275478e-09

# return analysis model
p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, return_analysis=TRUE)
#> lavaan 0.7-3 ended normally after 1 iteration
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         5
#> 
#>   Number of observations                            50
#> 
#> Model Test User Model:
#>                                                       
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0
#>                                                       
#>   Browne's residual (NT model-based) test             
#>   Test statistic                                 0.000
#>   Degrees of freedom                                 0

# data generation properties
N <- 1000
dat <- gen_mediation(n = N, a = .8, b = -.7, cprime = .2,
           sd.X = 2, sd.Y = 3, sd.M = 2)
descript(dat) # specific SDs
#> # A tibble: 3 × 11
#>   VARIABLE     n    mean    trim  median    sd   IQR     skew    kurt   min
#>   <fct>    <dbl>   <dbl>   <dbl>   <dbl> <dbl> <dbl>    <dbl>   <dbl> <dbl>
#> 1 X         1000  0.0381  0.0558  0.134   2.07  2.87 -0.0882  -0.112  -6.51
#> 2 M         1000  0.0441  0.0334  0.0837  2.05  2.77  0.00947  0.0319 -7.01
#> 3 Y         1000 -0.0906 -0.0790 -0.101   3.07  4.07 -0.0355  -0.156  -9.34
#> # ℹ 1 more variable: max <dbl>

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#> (Intercept)           X 
#>  0.01342498  0.80392087 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#> (Intercept)           M           X 
#> -0.06648038 -0.64186050  0.10933033 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#> -0.07509734 -0.40667472 

# same properties, but dichotomous X variable
dat <- gen_mediation(n = N, a = .8, b = -.7, cprime = .2,
           sd.X = 2, sd.Y = 3, sd.M = 2, dichotomous.X = TRUE)
descript(dat) # specific SDs
#> # A tibble: 3 × 11
#>   VARIABLE     n   mean   trim median    sd   IQR   skew    kurt    min   max
#>   <fct>    <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>  <dbl>   <dbl>  <dbl> <dbl>
#> 1 X         2000  2      2      2      2.00  4    0      -2.00     0     4   
#> 2 M         2000  1.56   1.55   1.51   2.02  3.25 0.0356 -0.859   -3.65  7.29
#> 3 Y         2000 -0.636 -0.642 -0.666  3.02  4.03 0.0402  0.0427 -11.1  11.0 

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#> (Intercept)           X 
#> -0.05499232  0.80602136 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#> (Intercept)           M           X 
#>   0.1352354  -0.6493776   0.1198056 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#>   0.1709461  -0.4036065 

# \donttest{

  # power to detect mediation
  p_mediation(n=50, a=sqrt(.35), b=sqrt(.35), cprime=.39) |>
    Spower(parallel=TRUE, replications=1000)
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 6
#>       n       a       b cprime sig.level power
#>   <dbl>   <dbl>   <dbl>  <dbl>     <dbl> <lgl>
#> 1    50 0.59161 0.59161   0.39      0.05 NA   
#> 
#> Estimate of power: 0.996
#> 95% Confidence Interval: [0.992, 1.000]
#> Execution time (H:M:S): 00:00:18

  # sample size estimate for .95 power
  p_mediation(n=interval(50,200), a=sqrt(.35), b=sqrt(.35), cprime=.39) |>
    Spower(power=.95, parallel=TRUE)
#> 
#> ── Spower Results ──────────────────────────────────────────────────────────────
#> 
#> Design conditions:
#> 
#> # A tibble: 1 × 6
#>       n       a       b cprime sig.level power
#>   <dbl>   <dbl>   <dbl>  <dbl>     <dbl> <dbl>
#> 1    NA 0.59161 0.59161   0.39      0.05  0.95
#> 
#> Estimate of n: 51.0
#> 95% Confidence Interval: [NA, NA]
#> Execution time (H:M:S): 00:19:51

# }
```

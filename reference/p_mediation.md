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

- two.tailed:

  logical; should a two-tailed or one-tailed test be used?

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
#> [1] 0.0008857608
p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, dichotomous.X=TRUE)
#> [1] 2.136133e-10

# return analysis model
p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, return_analysis=TRUE)
#> lavaan 0.6-21 ended normally after 1 iteration
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

# data generation properties
N <- 1000
dat <- gen_mediation(n = N, a = .8, b = -.7, cprime = .2,
           sd.X = 2, sd.Y = 3, sd.M = 2)
descript(dat) # specific SDs
#> # A tibble: 3 × 12
#>   VARS      n     mean    trim    sd   skew     kurt   min   P25     P50   P75
#>   <fct> <dbl>    <dbl>   <dbl> <dbl>  <dbl>    <dbl> <dbl> <dbl>   <dbl> <dbl>
#> 1 X      1000 -0.00310 -0.0220  2.04 0.134  -0.0358  -5.91 -1.42 -0.0496  1.42
#> 2 M      1000 -0.0198  -0.0408  2.08 0.118   0.00681 -6.07 -1.43 -0.0235  1.35
#> 3 Y      1000 -0.0975  -0.0941  2.99 0.0358  0.182   -9.01 -2.11 -0.103   1.91
#> # ℹ 1 more variable: max <dbl>

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#> (Intercept)           X 
#> -0.01725518  0.81214924 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#> (Intercept)           M           X 
#>  -0.1107392  -0.6987668   0.1870986 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#> -0.09868181 -0.38040435 

# same properties, but dichotomous X variable
dat <- gen_mediation(n = N, a = .8, b = -.7, cprime = .2,
           sd.X = 2, sd.Y = 3, sd.M = 2, dichotomous.X = TRUE)
descript(dat) # specific SDs
#> # A tibble: 3 × 12
#>   VARS      n   mean   trim    sd    skew    kurt    min     P25    P50   P75
#>   <fct> <dbl>  <dbl>  <dbl> <dbl>   <dbl>   <dbl>  <dbl>   <dbl>  <dbl> <dbl>
#> 1 X      2000  2      2      2.00  0      -2.00     0     0       2      4   
#> 2 M      2000  1.61   1.61   2.02 -0.0236 -0.854   -3.94 -0.0420  1.63   3.26
#> 3 Y      2000 -0.773 -0.779  3.01  0.0426  0.0297 -10.4  -2.84   -0.809  1.29
#> # ℹ 1 more variable: max <dbl>

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#> (Intercept)           X 
#> -0.02681293  0.81796198 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#> (Intercept)           M           X 
#>  0.01946498 -0.71296346  0.17727459 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#>  0.03858161 -0.40590242 

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
#> Estimate of power: 0.992
#> 95% Confidence Interval: [0.986, 0.998]
#> Execution time (H:M:S): 00:00:16

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
#> Estimate of n: 50.0
#> 95% Confidence Interval: [NA, NA]
#> Execution time (H:M:S): 00:18:16

# }
```

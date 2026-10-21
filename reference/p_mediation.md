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
#> [1] 9.312432e-05
p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, dichotomous.X=TRUE)
#> [1] 1.606742e-09

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
#>   VARIABLE     n    mean    trim   median    sd   IQR   skew    kurt   min   max
#>   <fct>    <dbl>   <dbl>   <dbl>    <dbl> <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl>
#> 1 X         1000 0.0132  0.00657 -0.00156  1.93  2.57 0.0620  0.188  -6.34  8.29
#> 2 M         1000 0.00402 0.00225  0.0381   1.95  2.61 0.0116 -0.0202 -6.32  6.11
#> 3 Y         1000 0.0863  0.0561  -0.0212   2.91  3.77 0.104   0.0626 -9.58  8.89

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#>  (Intercept)            X 
#> -0.006304465  0.783594684 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#> (Intercept)           M           X 
#>  0.08614592 -0.71280752  0.22940993 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#>  0.09063979 -0.32914225 

# same properties, but dichotomous X variable
dat <- gen_mediation(n = N, a = .8, b = -.7, cprime = .2,
           sd.X = 2, sd.Y = 3, sd.M = 2, dichotomous.X = TRUE)
descript(dat) # specific SDs
#> # A tibble: 3 × 11
#>   VARIABLE     n   mean   trim median    sd   IQR    skew   kurt    min   max
#>   <fct>    <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>   <dbl>  <dbl>  <dbl> <dbl>
#> 1 X         2000  2      2      2      2.00  4     0      -2.00    0     4   
#> 2 M         2000  1.63   1.64   1.65   2.00  3.24 -0.0140 -0.861  -3.75  6.43
#> 3 Y         2000 -0.761 -0.776 -0.716  2.92  3.75  0.0379  0.272 -12.9   9.96

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#> (Intercept)           X 
#>  0.01237286  0.81127665 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#> (Intercept)           M           X 
#> -0.06629501 -0.66002780  0.19224020 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#> -0.07446144 -0.34322494 

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
#> Estimate of power: 0.995
#> 95% Confidence Interval: [0.991, 0.999]
#> Execution time (H:M:S): 00:00:13

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
#> Execution time (H:M:S): 00:16:40

# }
```

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
#> [1] 2.177183e-07
p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, dichotomous.X=TRUE)
#> [1] 8.061746e-08

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
#> # A tibble: 3 × 14
#>   VARS      n  miss    mean trimmed    sd   mad skewness kurtosis   min  Q_25
#>   <fct> <dbl> <dbl>   <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl> <dbl> <dbl>
#> 1 X      1000    NA  0.112   0.131   2.12  2.18  -0.0612  -0.182  -6.59 -1.34
#> 2 M      1000    NA  0.0822  0.0963  2.02  2.11  -0.0466  -0.163  -5.67 -1.30
#> 3 Y      1000    NA -0.147  -0.137   2.97  3.01  -0.0611  -0.0395 -9.53 -2.19
#> # ℹ 3 more variables: Q_50 <dbl>, Q_75 <dbl>, max <dbl>

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#>  (Intercept)            X 
#> -0.004358922  0.771736024 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#> (Intercept)           M           X 
#>  -0.1116111  -0.6904618   0.1903087 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#>  -0.1086015  -0.3425456 

# same properties, but dichotomous X variable
dat <- gen_mediation(n = N, a = .8, b = -.7, cprime = .2,
           sd.X = 2, sd.Y = 3, sd.M = 2, dichotomous.X = TRUE)
descript(dat) # specific SDs
#> # A tibble: 3 × 14
#>   VARS      n  miss   mean trimmed    sd   mad skewness kurtosis    min    Q_25
#>   <fct> <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl>   <dbl>
#> 1 X      2000    NA  2       2      2.00  2.97   0       -2.00     0     0     
#> 2 M      2000    NA  1.60    1.61   1.96  2.35  -0.0348  -0.724   -4.68  0.0205
#> 3 Y      2000    NA -0.746  -0.738  3.06  3.03  -0.0256   0.0621 -12.4  -2.81  
#> # ℹ 3 more variables: Q_50 <dbl>, Q_75 <dbl>, max <dbl>

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#> (Intercept)           X 
#>  0.04366232  0.78023708 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#>  (Intercept)            M            X 
#>  0.003028295 -0.752976048  0.229479310 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#> -0.02984838 -0.35802052 

# \donttest{

  # power to detect mediation
  p_mediation(n=50, a=sqrt(.35), b=sqrt(.35), cprime=.39) |>
    Spower(parallel=TRUE, replications=1000)
#> 
#> Execution time (H:M:S): 00:00:21
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n cprime sig.level power
#>   <dbl>  <dbl>     <dbl> <lgl>
#> 1    50   0.39      0.05 NA   
#> 
#> Estimate of power: 0.996
#> 95% Confidence Interval: [0.992, 1.000]

  # sample size estimate for .95 power
  p_mediation(n=interval(50,200), a=sqrt(.35), b=sqrt(.35), cprime=.39) |>
    Spower(power=.95, parallel=TRUE)
#> 
#> Execution time (H:M:S): 00:24:11
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n cprime sig.level power
#>   <dbl>  <dbl>     <dbl> <dbl>
#> 1    NA   0.39      0.05  0.95
#> 
#> Estimate of n: 51.0
#> 95% Predicted Confidence Interval: [NA, NA]

# }
```

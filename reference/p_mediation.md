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
#> [1] 2.359283e-05
p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, dichotomous.X=TRUE)
#> [1] 2.584064e-08

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
#> 1 X      1000    NA  0.0384  0.0358  1.89  1.89  -0.0207   -0.143 -5.64 -1.28
#> 2 M      1000    NA  0.0716  0.0828  1.97  2.00  -0.0815   -0.137 -6.54 -1.26
#> 3 Y      1000    NA -0.0384 -0.0831  3.03  3.20   0.117    -0.227 -9.84 -2.18
#> # ℹ 3 more variables: Q_50 <dbl>, Q_75 <dbl>, max <dbl>

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#> (Intercept)           X 
#>  0.03962081  0.83270952 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#>  (Intercept)            M            X 
#>  0.003421063 -0.750061890  0.308934171 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#>  -0.0262970  -0.3156495 

# same properties, but dichotomous X variable
dat <- gen_mediation(n = N, a = .8, b = -.7, cprime = .2,
           sd.X = 2, sd.Y = 3, sd.M = 2, dichotomous.X = TRUE)
descript(dat) # specific SDs
#> # A tibble: 3 × 14
#>   VARS      n  miss   mean trimmed    sd   mad skewness kurtosis    min   Q_25
#>   <fct> <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl>  <dbl>
#> 1 X      2000    NA  2       2      2.00  2.97   0        -2.00    0     0    
#> 2 M      2000    NA  1.56    1.56   2.03  2.47   0.0213   -0.903  -3.87 -0.124
#> 3 Y      2000    NA -0.764  -0.795  3.02  2.94   0.133     0.144 -10.5  -2.77 
#> # ℹ 3 more variables: Q_50 <dbl>, Q_75 <dbl>, max <dbl>

# two-step regression-based estimates (not used)
lm(M ~ X, data=dat) |> coef()       # a
#> (Intercept)           X 
#> -0.07499256  0.81981342 
lm(Y ~ M + X, data=dat) |> coef()   # b and cprime
#>  (Intercept)            M            X 
#> -0.003718004 -0.632247676  0.114442992 
lm(Y ~ X, data=dat) |> coef()       # c = cprime + a*b
#> (Intercept)           X 
#>  0.04369587 -0.40388213 

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
#> Estimate of power: 0.998
#> 95% Confidence Interval: [0.995, 1.000]

  # sample size estimate for .95 power
  p_mediation(n=interval(50,200), a=sqrt(.35), b=sqrt(.35), cprime=.39) |>
    Spower(power=.95, parallel=TRUE)
#> 
#> Execution time (H:M:S): 00:24:07
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n cprime sig.level power
#>   <dbl>  <dbl>     <dbl> <dbl>
#> 1    NA   0.39      0.05  0.95
#> 
#> Estimate of n: 56.3
#> 95% Predicted Confidence Interval: [NA, 51.0]

# }
```

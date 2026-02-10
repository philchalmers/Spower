# p-value from (generalized) linear regression model simulations with fixed predictors

p-values associated with (generalized) linear regression model. Requires
a pre-specified design matrix (`X`).

## Usage

``` r
p_glm(
  formula,
  X,
  betas,
  test,
  sigma = NULL,
  family = gaussian(),
  gen_fun = gen_glm,
  return_analysis = FALSE,
  ...
)

gen_glm(formula, X, betas, sigma = NULL, family = gaussian(), ...)
```

## Arguments

- formula:

  formula passed to either [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`glm`](https://rdrr.io/r/stats/glm.html)

- X:

  a data.frame containing the covariates

- betas:

  vector of slope coefficients that match the `model.matrix` version of
  `X`

- test:

  character vector specifying the test to pass to
  [`lht`](https://rdrr.io/pkg/car/man/linearHypothesis.html). Can also
  be a list of character vectors to evaluate multiple tests

- sigma:

  residual standard deviation for linear model. Only used when
  `family = 'gaussian'`

- family:

  family of distributions to use (see
  [`family`](https://rdrr.io/r/stats/family.html))

- gen_fun:

  function used to generate the required discrete data. Object returned
  must be a `data.frame`. Default uses `gen_glm`. User defined version
  of this function must include the argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

## Value

a single p-value

## See also

[`p_lm.R2`](https://philchalmers.github.io/Spower/reference/p_lm.R2.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
X <- data.frame(G = factor(rep(c('control', 'treatment'), each=50)),
                C = sample(50:100, 100, replace=TRUE))
head(X)
#>         G  C
#> 1 control 96
#> 2 control 72
#> 3 control 83
#> 4 control 66
#> 5 control 54
#> 6 control 85

# ANCOVA setup
p_glm(y ~ G + C, test="Gtreatment = 0",
  X=X, betas=c(10, .3, 1), sigma=1)
#> [1] 0.02103663

# return analysis model
p_glm(y ~ G + C, test="Gtreatment = 0",
  X=X, betas=c(10, .3, 1), sigma=1, return_analysis=TRUE)
#> 
#> Call:
#> lm(formula = formula, data = X)
#> 
#> Coefficients:
#> (Intercept)   Gtreatment            C  
#>      9.7051       0.5201       1.0033  
#> 


# ANCOVA setup with logistic regression
p_glm(y ~ G + C, test="Gtreatment = 0",
  X=X, betas=c(-2, .5, .01), family=binomial())
#> [1] 0.4833122

# ANCOVA setup with poisson regression
p_glm(y ~ G + C, test="Gtreatment = 0",
  X=X, betas=c(-2, .5, .01), family=poisson())
#> [1] 0.01291868

# \donttest{

# test whether two slopes differ given different samples.
#   To do this setup data as an MLR where a binary variable S
#   is used to reflect the second sample, and the interaction
#   effect evaluates the magnitude of the slope difference
gen_twogroup <- function(n, dbeta, sdx1, sdx2, sigma, n2_n1 = 1, ...){
  X1 <- rnorm(n, sd=sdx1)
  X2 <- rnorm(n*n2_n1, sd=sdx2)
  X <- c(X1, X2)
  N <- length(X)
  S <- c(rep(0, n), rep(1, N-n))
  y <- dbeta * X*S + rnorm(N, sd=sigma)
  dat <- data.frame(y, X, S)
  dat
}

# prospective power using test that interaction effect is equal to 0
p_glm(formula=y~X*S, test="X:S = 0",
    n=100, sdx1=1, sdx2=2, dbeta=0.2,
    sigma=0.5, gen_fun=gen_twogroup) |> Spower(replications=1000)
#> 
#> Execution time (H:M:S): 00:00:02
#> Design conditions: 
#> 
#> # A tibble: 1 × 8
#>   test    sigma     n  sdx1  sdx2 dbeta sig.level power
#>   <chr>   <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1 X:S = 0   0.5   100     1     2   0.2      0.05 NA   
#> 
#> Estimate of power: 0.923
#> 95% Confidence Interval: [0.906, 0.940]

# }

```

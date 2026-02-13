# p-value from global linear regression model simulation

p-values associated with linear regression model using fixed/random
independent variables. Focus is on the omnibus behavior of the R^2
statistic.

## Usage

``` r
p_lm.R2(
  n,
  R2,
  k,
  R2_0 = 0,
  k.R2_0 = 0,
  R2.resid = 1 - R2,
  fixed = TRUE,
  return_analysis = FALSE,
  ...
)
```

## Arguments

- n:

  sample size

- R2:

  R-squared effect size

- k:

  number of IVs

- R2_0:

  null hypothesis for R-squared

- k.R2_0:

  number of IVs associated with the null hypothesis model

- R2.resid:

  residual R-squared value, typically used when comparing nested models
  when fit sequentially (e.g., comparing model A vs B when model
  involves the structure A -\> B -\> C)

- fixed:

  logical; if FALSE then the data are random generated according to a
  joint multivariate normal distribution

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

## Value

a single p-value

## See also

[`p_glm`](https://philchalmers.github.io/Spower/reference/p_glm.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# 5 fixed IVs, R^2 = .1, sample size of 95
p_lm.R2(n=95, R2=.1, k=5)
#> [1] 0.004587803

# random model
p_lm.R2(n=95, R2=.1, k=5, fixed=FALSE)
#> [1] 0.001075562

# return analysis model
p_lm.R2(n=95, R2=.1, k=5, return_analysis=TRUE)
#> 
#> Call:
#> lm(formula = y ~ ., data = df)
#> 
#> Coefficients:
#> (Intercept)           X1           X2           X3           X4           X5  
#>     0.05595      0.42595      0.07513      0.08230     -0.05399      0.15285  
#> 
```

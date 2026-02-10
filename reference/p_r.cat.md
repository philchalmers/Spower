# p-value from tetrachoric/polychoric or polyserial

Generates correlated X-Y data and returns a p-value to assess the null
of no correlation in the population. The X-Y data are generated assuming
a multivariate normal distribution and subsequently discretized for one
or both of the variables.

## Usage

``` r
p_r.cat(
  n,
  r,
  tauX,
  rho = 0,
  tauY = NULL,
  ML = TRUE,
  two.tailed = TRUE,
  score = FALSE,
  gen_fun = gen_r,
  return_analysis = FALSE,
  ...
)
```

## Arguments

- n:

  sample size

- r:

  correlation prior to the discretization (recovered via the
  polyserial/polychoric estimates)

- tauX:

  intercept parameters used for discretizing the X variable

- rho:

  population coefficient to test against

- tauY:

  intercept parameters used for discretizing the Y variable. If missing
  a polyserial correlation will be estimated, otherwise a
  tetrachoric/polychoric correlation will be estimated

- ML:

  logical; use maximum-likelihood estimation?

- two.tailed:

  logical; should a two-tailed or one-tailed test be used?

- score:

  logical; should the SE be based at the null hypothesis (score test) or
  the ML estimate (Wald test)? The former is the canonical form for a
  priori power analyses though requires twice as many computations as
  the Wald test approach

- gen_fun:

  function used to generate the required continuous bivariate data
  (prior to truncation). Object returned must be a `matrix` with two
  columns. Default uses
  [`gen_r`](https://philchalmers.github.io/Spower/reference/p_r.md) to
  generate conditionally dependent data from a bivariate normal
  distribution. User defined version of this function must include the
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

[`gen_r`](https://philchalmers.github.io/Spower/reference/p_r.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# 100 observations, .5 correlation, tetrachoric estimate
p_r.cat(100, r=.5, tauX=0, tauY=1)
#> [1] 0.02083308

# return analysis model
p_r.cat(100, r=.5, tauX=0, tauY=1, return_analysis=TRUE)
#> 
#> Polychoric Correlation, ML est. = 0.4779 (0.1793)
#> 
#>   Row Threshold
#>   Threshold Std.Err.
#>       1.175   0.1624
#> 
#> 
#>   Column Threshold
#>   Threshold Std.Err.
#>    -0.07527   0.1255

# Wald test
p_r.cat(100, r=.5, tauX=0, tauY=1, score=FALSE)
#> [1] 0.01119212

# polyserial estimate (Y continuous)
p_r.cat(50, r=.5, tauX=0)
#> [1] 0.006380214
```

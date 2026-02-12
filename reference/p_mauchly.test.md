# p-value from Mauchly's Test of Sphericity simulation

Perform simulation experiment for Mauchly's Test of Sphericity using the
function `mauchlys.test`, returning a p-value. Assumes the data are from
a multivariate normal distribution, however this can be modified.

## Usage

``` r
p_mauchly.test(
  n,
  sigma,
  gen_fun = gen_mauchly.test,
  return_analysis = FALSE,
  ...
)

gen_mauchly.test(n, sigma, ...)

mauchlys.test(X)
```

## Arguments

- n:

  sample size

- sigma:

  symmetric covariance/correlation matrix passed to `gen_fun`

- gen_fun:

  function used to generate the required data. Object returned must be a
  `matrix` with `K` columns and `n` rows. Default uses
  `gen_mauchly.test` to generate multivariate normal samples. User
  defined version of this function must include the argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

- X:

  a matrix with `k` columns and `n` rows

## Value

a single p-value

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
sigma <- diag(c(1,2,1))
sigma
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    2    0
#> [3,]    0    0    1

# H0 test that sphericity holds
p_mauchly.test(100, sigma=sigma)
#> [1] 0.06869854

# return analysis object
p_mauchly.test(100, sigma=sigma, return_analysis=TRUE)
#>           W df    p.value
#> 1 0.9523496  2 0.09141641

# Null is true
sigma.H0 <- diag(3)
p_mauchly.test(100, sigma=sigma.H0)
#> [1] 0.04578572


# \donttest{
    # empirical power estimate
    p_mauchly.test(100, sigma=sigma) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:09
#> Design conditions: 
#> 
#> # A tibble: 1 × 3
#>       n sig.level power
#>   <dbl>     <dbl> <lgl>
#> 1   100      0.05 NA   
#> 
#> Estimate of power: 0.598
#> 95% Confidence Interval: [0.589, 0.608]

    # empirical Type I error estimate
    p_mauchly.test(100, sigma=sigma.H0) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:09
#> Design conditions: 
#> 
#> # A tibble: 1 × 3
#>       n sig.level power
#>   <dbl>     <dbl> <lgl>
#> 1   100      0.05 NA   
#> 
#> Estimate of power: 0.052
#> 95% Confidence Interval: [0.047, 0.056]
# }
```

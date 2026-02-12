# p-value from proportion test simulation

Generates single and multi-sample data for proportion tests and return a
p-value. Uses [`binom.test`](https://rdrr.io/r/stats/binom.test.html)
for one-sample applications and
[`prop.test`](https://rdrr.io/r/stats/prop.test.html) otherwise.

## Usage

``` r
p_prop.test(
  n,
  h,
  prop = NULL,
  pi = 0.5,
  n.ratios = rep(1, length(prop)),
  two.tailed = TRUE,
  correct = TRUE,
  exact = FALSE,
  gen_fun = gen_prop.test,
  return_analysis = FALSE,
  ...
)

gen_prop.test(
  n,
  h,
  prop = NULL,
  pi = 0.5,
  n.ratios = rep(1, length(prop)),
  ...
)
```

## Arguments

- n:

  sample size per group

- h:

  Cohen's h effect size; only supported for one-sample analysis.

  Note that it's important to specify the null value `pi` when supplying
  this effect size as the power changes depending on these specific
  values (see example below).

- prop:

  sample probability/proportions of success. If a vector with two-values
  or more elements are supplied then a multi-samples test will be used.
  Matrices are also supported

- pi:

  probability of success to test against (default is .5). Ignored for
  two-sample tests

- n.ratios:

  allocation ratios reflecting the sample size ratios. Default of 1 sets
  the groups to be the same size (n \* n.ratio)

- two.tailed:

  logical; should a two-tailed or one-tailed test be used?

- correct:

  logical; use Yates' continuity correction?

- exact:

  logical; use fisher's exact test via
  [`fisher.test`](https://rdrr.io/r/stats/fisher.test.html)? Use of this
  flag requires that `prop` was specified as a matrix

- gen_fun:

  function used to generate the required discrete data. Object returned
  must be a `matrix` with two rows and 1 or more columns. Default uses
  `gen_prop.test`. User defined version of this function must include
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

`gen_prop.test`

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# one sample, 50 observations, tested against pi = .5 by default
p_prop.test(50, prop=.65)
#> [1] 0.03283914

# return analysis model
p_prop.test(50, prop=.65, return_analysis = TRUE)
#> 
#>  Exact binomial test
#> 
#> data:  dat[1, 1] and n
#> number of successes = 33, number of trials = 50, p-value = 0.03284
#> alternative hypothesis: true probability of success is not equal to 0.5
#> 95 percent confidence interval:
#>  0.5123475 0.7879453
#> sample estimates:
#> probability of success 
#>                   0.66 
#> 

# specified using h and pi
h <- pwr::ES.h(.65, .4)
p_prop.test(50, h=h, pi=.4)
#> [1] 0.05939392
p_prop.test(50, h=-h, pi=.65)
#> [1] 6.656485e-07

# two-sample test
p_prop.test(50, prop=c(.5, .65))
#> [1] 0.04511188

# two-sample test, unequal ns
p_prop.test(50, prop=c(.5, .65), n.ratios = c(1,2))
#> [1] 0.03167776

# three-sample test, group2 twice as large as others
p_prop.test(50, prop=c(.5, .65, .7), n.ratios=c(1,2,1))
#> [1] 0.02882835

# Fisher exact test
p_prop.test(50, prop=matrix(c(.5, .65, .7, .5), 2, 2))
#> [1] 0.03072934

# \donttest{
    # compare simulated results to pwr package

    # one-sample tests
    (h <- pwr::ES.h(0.5, 0.4))
#> [1] 0.2013579
    pwr::pwr.p.test(h=h, n=60)
#> 
#>      proportion power calculation for binomial distribution (arcsine transformation) 
#> 
#>               h = 0.2013579
#>               n = 60
#>       sig.level = 0.05
#>           power = 0.3447014
#>     alternative = two.sided
#> 

    # uses binom.test (need to specify null location as this matters!)
    Spower(p_prop.test(n=60, h=h, pi=.4))
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n    pi sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    60   0.4      0.05 NA   
#> 
#> Estimate of power: 0.354
#> 95% Confidence Interval: [0.345, 0.363]
    Spower(p_prop.test(n=60, prop=.5, pi=.4))
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n  prop    pi sig.level power
#>   <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1    60   0.5   0.4      0.05 NA   
#> 
#> Estimate of power: 0.358
#> 95% Confidence Interval: [0.348, 0.367]

    # compare with switched null
    Spower(p_prop.test(n=60, h=h, pi=.5))
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n    pi sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    60   0.5      0.05 NA   
#> 
#> Estimate of power: 0.266
#> 95% Confidence Interval: [0.257, 0.274]
    Spower(p_prop.test(n=60, prop=.4, pi=.5))
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n  prop    pi sig.level power
#>   <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1    60   0.4   0.5      0.05 NA   
#> 
#> Estimate of power: 0.248
#> 95% Confidence Interval: [0.239, 0.256]

    # two-sample test, one-tailed
    (h <- pwr::ES.h(0.67, 0.5))
#> [1] 0.3469169
    pwr::pwr.2p.test(h=h, n=80, alternative="greater")
#> 
#>      Difference of proportion power calculation for binomial distribution (arcsine transformation) 
#> 
#>               h = 0.3469169
#>               n = 80
#>       sig.level = 0.05
#>           power = 0.7085801
#>     alternative = greater
#> 
#> NOTE: same sample sizes
#> 
    p_prop.test(n=80, prop=c(.67, .5), two.tailed=FALSE,
      correct=FALSE) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n two.tailed correct sig.level power
#>   <dbl> <lgl>      <lgl>       <dbl> <lgl>
#> 1    80 FALSE      FALSE        0.05 NA   
#> 
#> Estimate of power: 0.698
#> 95% Confidence Interval: [0.689, 0.707]

    # same as above, but with continuity correction (default)
    p_prop.test(n=80, prop=c(.67, .5), two.tailed=FALSE) |>
      Spower()
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n two.tailed sig.level power
#>   <dbl> <lgl>          <dbl> <lgl>
#> 1    80 FALSE           0.05 NA   
#> 
#> Estimate of power: 0.634
#> 95% Confidence Interval: [0.625, 0.644]

    # three-sample joint test, equal n's
    p_prop.test(n=50, prop=c(.6,.4,.7)) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:03
#> Design conditions: 
#> 
#> # A tibble: 1 × 3
#>       n sig.level power
#>   <dbl>     <dbl> <lgl>
#> 1    50      0.05 NA   
#> 
#> Estimate of power: 0.806
#> 95% Confidence Interval: [0.798, 0.813]

# }
```

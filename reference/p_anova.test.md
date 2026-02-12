# p-value from one-way ANOVA simulation

Generates continuous multi-sample data to be analyzed by a one-way
ANOVA, and return a p-value. Uses the function
[`oneway.test`](https://rdrr.io/r/stats/oneway.test.html) to perform the
analyses. The data and associated test assume that the conditional
observations are normally distributed and have have equal variance by
default, however these may be modified.

## Usage

``` r
p_anova.test(
  n,
  k,
  f,
  n.ratios = rep(1, k),
  two.tailed = TRUE,
  var.equal = TRUE,
  means = NULL,
  sds = NULL,
  gen_fun = gen_anova.test,
  return_analysis = FALSE,
  ...
)

gen_anova.test(n, k, f, n.ratios = rep(1, k), means = NULL, sds = NULL, ...)
```

## Arguments

- n:

  sample size per group

- k:

  number of groups

- f:

  Cohen's f effect size

- n.ratios:

  allocation ratios reflecting the sample size ratios. Default of 1 sets
  the groups to be the same size (n \* n.ratio)

- two.tailed:

  logical; should a two-tailed or one-tailed test be used?

- var.equal:

  logical; use the pooled SE estimate instead of the Welch correction
  for unequal variances?

- means:

  (optional) vector of means. When specified the input `f` is ignored

- sds:

  (optional) vector of SDs. When specified the input `f` is ignored

- gen_fun:

  function used to generate the required data. Object returned must be a
  `matrix` with k rows and k columns of numeric data. Default uses
  `gen_anova.test`. User defined version of this function must include
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

`gen_anova.test`

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# n=50 in 3 groups, "medium" effect size
p_anova.test(50, k=3, f=.25)
#> [1] 0.0008046224

# return analysis model
p_anova.test(50, k=3, f=.25, return_analysis=TRUE)
#> 
#>  One-way analysis of means
#> 
#> data:  DV and group
#> F = 4.2735, num df = 2, denom df = 147, p-value = 0.0157
#> 

# explicit means/sds
p_anova.test(50, 3, means=c(0,0,1), sds=c(1,2,1))
#> [1] 3.392539e-05

# \donttest{
  # compare simulated results to pwr package
  pwr::pwr.anova.test(f=0.28, k=4, n=20)
#> 
#>      Balanced one-way analysis of variance power calculation 
#> 
#>               k = 4
#>               n = 20
#>               f = 0.28
#>       sig.level = 0.05
#>           power = 0.5149793
#> 
#> NOTE: n is number in each group
#> 
  p_anova.test(n=20, k=4, f=.28) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:09
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n     k     f sig.level power
#>   <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1    20     4  0.28      0.05 NA   
#> 
#> Estimate of power: 0.511
#> 95% Confidence Interval: [0.502, 0.521]
# }
```

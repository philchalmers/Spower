# p-value from Wilcoxon (signed rank) test simulation

Simulates data given one (Wilcoxon) or two (Mann-Whitney) parent
distributions and returns a p-value. Can also be used for power analyses
related to sign tests.

## Usage

``` r
p_wilcox.test(
  n,
  d,
  n2_n1 = 1,
  mu = 0,
  type = c("two.sample", "one.sample", "paired"),
  exact = NULL,
  correct = TRUE,
  two.tailed = TRUE,
  parent1 = function(n, d) rnorm(n, d, 1),
  parent2 = function(n, d) rnorm(n, 0, 1),
  return_analysis = FALSE
)
```

## Arguments

- n:

  sample size per group. For paired samples this corresponds to the
  number of pairings (hence, half of the data points observed)

- d:

  effect size passed to `parent` functions

- n2_n1:

  sample size ratio

- mu:

  parameter used to form the null hypothesis

- type:

  type of analysis to use (two-sample, one-sample, or paired)

- exact:

  a logical indicating whether an exact p-value should be computed

- correct:

  a logical indicating whether to apply continuity correction in the
  normal approximation for the p-value

- two.tailed:

  logical; use two-tailed test?

- parent1:

  data generation function for first group. Ideally should have SDs = 1
  so that `d` reflects a standardized difference

- parent2:

  same as `parent1`, but for the second group

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

## Value

a single p-value

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# with normal distributions defaults d is standardized
p_wilcox.test(100, .5)
#> [1] 0.01981824
p_wilcox.test(100, .5, type = 'paired')  # n = number of pairs
#> [1] 0.0006436109
p_wilcox.test(100, .5, type = 'one.sample')
#> [1] 2.782262e-08

# return analysis object
p_wilcox.test(100, .5, return_analysis = TRUE)
#> 
#>  Wilcoxon rank sum test with continuity correction
#> 
#> data:  dat1 and dat2
#> W = 6154, p-value = 0.004826
#> alternative hypothesis: true location shift is not equal to 0
#> 

# using chi-squared distributions (standardizing to 0-1)
p_wilcox.test(100, .5, type = 'one.sample',
   parent1 = function(n, d) rchisq(n, df=10) - 10 + d)
#> [1] 0.6486934
p_wilcox.test(100, .5,
   parent1 = function(n, d) (rchisq(n, df=10) - 10)/sqrt(20) + d,
   parent2 = function(n, d) (rchisq(n, df=10) - 10)/sqrt(20))
#> [1] 0.01844069
```

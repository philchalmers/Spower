# p-value from Kruskal-Wallis Rank Sum Test simulation

Simulates data given two or more parent distributions and returns a
p-value using
[`kruskal.test`](https://rdrr.io/r/stats/kruskal.test.html). Default
generates data from Gaussian distributions, however this can be
modified.

## Usage

``` r
p_kruskal.test(
  n,
  k,
  means,
  n.ratios = rep(1, k),
  gen_fun = gen_kruskal.test,
  return_analysis = FALSE,
  ...
)

gen_kruskal.test(n, k, n.ratios, means, ...)
```

## Arguments

- n:

  sample size per group

- k:

  number of groups

- means:

  vector of means to control location parameters

- n.ratios:

  allocation ratios reflecting the sample size ratios. Default of 1 sets
  the groups to be the same size (n \* n.ratio)

- gen_fun:

  function used to generate the required data. Object returned must be a
  `list` of length `k`, where each element contains the sample data in
  each group. Default uses `gen_kruskal.test`. User defined version of
  this function must include the argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to pass to `gen_fun`

## Value

a single p-value

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# three group test where data generate from Gaussian distributions
p_kruskal.test(n=30, k=3, means=c(0, .5, .6))
#> [1] 0.004637737

# return analysis model
p_kruskal.test(n=30, k=3, means=c(0, .5, .6), return_analysis=TRUE)
#> 
#>  Kruskal-Wallis rank sum test
#> 
#> data:  dat
#> Kruskal-Wallis chi-squared = 4.5575, df = 2, p-value = 0.1024
#> 

# generate data from chi-squared distributions with different variances
gen_chisq <- function(n, k, n.ratios, means, dfs, ...){
  dat <- vector('list', k)
  ns <- n * n.ratios
  for(g in 1:k)
   dat[[g]] <- rchisq(ns[g], df=dfs[g]) - dfs[g] + means[g]
  dat
}

p_kruskal.test(n=30, k=3, means=c(0, 1, 2),
   gen_fun=gen_chisq, dfs=c(10, 15, 20))
#> [1] 0.1179166

# \donttest{
  # empirical power estimate
  p_kruskal.test(n=30, k=3, means=c(0, .5, .6)) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:06
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n     k sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    30     3      0.05 NA   
#> 
#> Estimate of power: 0.552
#> 95% Confidence Interval: [0.542, 0.561]
  p_kruskal.test(n=30, k=3, means=c(0, 1, 2), gen_fun=gen_chisq,
         dfs = c(10, 15, 20)) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:07
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n     k sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1    30     3      0.05 NA   
#> 
#> Estimate of power: 0.194
#> 95% Confidence Interval: [0.186, 0.202]

# }
```

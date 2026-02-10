# p-value from McNemar test simulation

Generates two-dimensional sample data for McNemar test and return a
p-value. Uses
[`mcnemar.test`](https://rdrr.io/r/stats/mcnemar.test.html).

## Usage

``` r
p_mcnemar.test(
  n,
  prop,
  OR = NULL,
  prop.disc = NULL,
  two.tailed = TRUE,
  correct = TRUE,
  gen_fun = gen_mcnemar.test,
  return_analysis = FALSE,
  ...
)

gen_mcnemar.test(n, prop, ...)
```

## Arguments

- n:

  total sample size

- prop:

  two-dimensional matrix of proportions/probabilities

- OR:

  instead of supplying the `prop` table, the odds ratio can be specified
  instead \\\pi\_{12}/\pi\_{21}\\. Also requires proportion of
  discordant pairings to be specified

- prop.disc:

  proportion of discordant pairings, \\\pi\_{12} + \pi\_{21}\\

- two.tailed:

  logical; should a two-tailed or one-tailed test be used?

- correct:

  logical; use continuity correction? Only applicable for 2x2 tables

- gen_fun:

  function used to generate the required discrete data. Object returned
  must be a `matrix` with k rows and k columns of counts. Default uses
  `gen_mcnemar.test`. User defined version of this function must include
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

`gen_mcnemar.test`

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# from ?mcnemar.test
Performance <- matrix(c(794, 86, 150, 570),
       nrow = 2,
       dimnames = list("1st Survey" = c("Approve", "Disapprove"),
                   "2nd Survey" = c("Approve", "Disapprove")))
(prop <- prop.table(Performance))
#>             2nd Survey
#> 1st Survey   Approve Disapprove
#>   Approve    0.49625    0.09375
#>   Disapprove 0.05375    0.35625

# one sample + test and resulting p-value
p_mcnemar.test(n=sum(Performance), prop=prop)
#> [1] 9.682145e-11

# return analysis model
p_mcnemar.test(n=sum(Performance), prop=prop, return_analysis=TRUE)
#> 
#>  McNemar's Chi-squared test with continuity correction
#> 
#> data:  dat
#> McNemar's chi-squared = 16.676, df = 1, p-value = 4.433e-05
#> 

# \donttest{

# post-hoc power (not recommended)
Spower(p_mcnemar.test(n=sum(Performance), prop=prop))
#> Error in eval(x, parent.frame()): object 'Performance' not found

# odds ratio + discordant proportions supplied instead
OR <- prop[1,2] / prop[2,1]
disc <- prop[1,2] + prop[2,1]
p_mcnemar.test(n=50, OR=.25, prop.disc=disc, two.tailed=FALSE) |>
  Spower(replications=30000)
#> 
#> Execution time (H:M:S): 00:00:05
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n    OR two.tailed sig.level power
#>   <dbl> <dbl> <lgl>          <dbl> <lgl>
#> 1    50  0.25 FALSE           0.05 NA   
#> 
#> Estimate of power: 0.337
#> 95% Confidence Interval: [0.332, 0.342]

# }
```

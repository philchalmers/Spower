# p-value from independent/paired samples t-test simulation

Generates one or two sets of continuous data group-level data according
to Cohen's effect size 'd', and returns a p-value. The data and
associated t-test assume that the conditional observations are normally
distributed and have have equal variance by default, however these may
be modified.

## Usage

``` r
p_t.test(
  n,
  d,
  mu = 0,
  r = NULL,
  type = "two.sample",
  n2_n1 = 1,
  two.tailed = TRUE,
  var.equal = TRUE,
  means = NULL,
  sds = NULL,
  conf.level = 0.95,
  gen_fun = gen_t.test,
  return_analysis = FALSE,
  ...
)

gen_t.test(
  n,
  d,
  n2_n1 = 1,
  r = NULL,
  type = "two.sample",
  means = NULL,
  sds = NULL,
  ...
)
```

## Arguments

- n:

  sample size per group, assumed equal across groups. For paired samples
  this corresponds to the number of pairs (hence, half the number of
  data points observed)

- d:

  Cohen's standardized effect size `d`. For the generated data this
  standardized mean appears in the first group (two-sample)/first time
  point (paired samples)

- mu:

  population mean to test against

- r:

  (optional) instead of specifying `d` specify a point-biserial
  correlation. Internally this is transformed into a suitable `d` value
  for the power computations

- type:

  type of t-test to use; can be `'two.sample'`, `'one.sample'`, or
  `'paired'`

- n2_n1:

  allocation ratio reflecting the same size ratio. Default of 1 sets the
  groups to be the same size. Only applicable when `type = 'two.sample'`

- two.tailed:

  logical; should a two-tailed or one-tailed test be used?

- var.equal:

  logical; use the classical or Welch corrected t-test?

- means:

  (optional) vector of means for each group. When specified the input
  `d` is ignored

- sds:

  (optional) vector of SDs for each group. If not specified and `d` is
  used then these are set to a vector of 1's

- conf.level:

  confidence interval level passed to
  [`t.test`](https://rdrr.io/r/stats/t.test.html)

- gen_fun:

  function used to generate the required two-sample data. Object
  returned must be a `list` containing one (one-sample) or two
  (independent samples/paired samples) elements, both of which are
  `numeric` vectors. Default uses `gen_t.test` to generate conditionally
  Gaussian distributed samples. User defined version of this function
  must include the argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

## Value

a single p-value

## See also

`gen_t.test`

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# sample size of 50 per group, "medium" effect size
p_t.test(n=50, d=0.5)
#> [1] 0.005589963

# point-biserial correlation effect size
p_t.test(n=50, r=.3)
#> [1] 0.004396717

# second group 2x as large as the first group
p_t.test(n=50, d=0.5, n2_n1 = 2)
#> [1] 0.03033605

# specify mean/SDs explicitly
p_t.test(n=50, means = c(0,1), sds = c(2,2))
#> [1] 0.009252317

# paired and one-sample tests
p_t.test(n=50, d=0.5, type = 'paired') # n = number of pairs
#> [1] 0.06232847
p_t.test(n=50, d=0.5, type = 'one.sample')
#> [1] 0.0007065882

# return analysis object
p_t.test(n=50, d=0.5, return_analysis=TRUE)
#> 
#>  Two Sample t-test
#> 
#> data:  dat[[1]] and dat[[2]]
#> t = 2.0981, df = 98, p-value = 0.03846
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  0.02490457 0.89441223
#> sample estimates:
#>  mean of x  mean of y 
#> 0.54009219 0.08043378 
#> 

# \donttest{
  # compare simulated results to pwr package

  pwr::pwr.t.test(d=0.2, n=60, sig.level=0.10,
             type="one.sample", alternative="two.sided")
#> 
#>      One-sample t test power calculation 
#> 
#>               n = 60
#>               d = 0.2
#>       sig.level = 0.1
#>           power = 0.4555818
#>     alternative = two.sided
#> 
  p_t.test(n=60, d=0.2, type = 'one.sample', two.tailed=TRUE) |>
         Spower(sig.level=.10)
#> 
#> Execution time (H:M:S): 00:00:02
#> Design conditions: 
#> 
#> # A tibble: 1 × 6
#>       n     d type       two.tailed sig.level power
#>   <dbl> <dbl> <chr>      <lgl>          <dbl> <lgl>
#> 1    60   0.2 one.sample TRUE             0.1 NA   
#> 
#> Estimate of power: 0.451
#> 95% Confidence Interval: [0.441, 0.460]

  pwr::pwr.t.test(d=0.3, power=0.80, type="two.sample",
                  alternative="greater")
#> 
#>      Two-sample t test power calculation 
#> 
#>               n = 138.0716
#>               d = 0.3
#>       sig.level = 0.05
#>           power = 0.8
#>     alternative = greater
#> 
#> NOTE: n is number in *each* group
#> 
  p_t.test(n=NA, d=0.3, type='two.sample', two.tailed=FALSE) |>
         Spower(power=0.80, interval=c(10,200))
#> 
#> Execution time (H:M:S): 00:00:20
#> Design conditions: 
#> 
#> # A tibble: 1 × 6
#>       n     d type       two.tailed sig.level power
#>   <dbl> <dbl> <chr>      <lgl>          <dbl> <dbl>
#> 1    NA   0.3 two.sample FALSE           0.05   0.8
#> 
#> Estimate of n: 136.9
#> 95% Predicted Confidence Interval: [135.2, 138.7]

# }


###### Custom data generation function

# Generate data such that:
#   - group 1 is from a negatively distribution (reversed X2(10)),
#   - group 2 is from a positively skewed distribution (X2(5))
#   - groups have equal variance, but differ by d = 0.5

args(gen_t.test)   ## can use these arguments as a basis, though must include ...
#> function (n, d, n2_n1 = 1, r = NULL, type = "two.sample", means = NULL, 
#>     sds = NULL, ...) 
#> NULL

# arguments df1 and df2 added; unused arguments caught within ...
my.gen_fun <- function(n, d, df1, df2, ...){
    group1 <- -1 * rchisq(n, df=df1)
       group2 <- rchisq(n, df=df2)
       # scale groups first given moments of the chi-square distribution,
       #   then add std mean difference
       group1 <- ((group1 + df1) / sqrt(2*df1))
       group2 <- ((group2 - df2) / sqrt(2*df2)) + d
       dat <- list(group1, group2)
       dat
}

# check the sample data properties
dat <- my.gen_fun(n=10000, d=.5, df1=10, df2=5)
sapply(dat, mean)
#> [1] -0.00159157  0.49555302
sapply(dat, sd)
#> [1] 0.9914658 0.9927854

p_t.test(n=100, d=0.5, gen_fun=my.gen_fun, df1=10, df2=5)
#> [1] 0.004709083

# \donttest{

  # power given Gaussian distributions
  p_t.test(n=100, d=0.5) |> Spower(replications=30000)
#> 
#> Execution time (H:M:S): 00:00:08
#> Design conditions: 
#> 
#> # A tibble: 1 × 4
#>       n     d sig.level power
#>   <dbl> <dbl>     <dbl> <lgl>
#> 1   100   0.5      0.05 NA   
#> 
#> Estimate of power: 0.942
#> 95% Confidence Interval: [0.939, 0.944]

  # estimate power given the customized data generating function
  p_t.test(n=100, d=0.5, gen_fun=my.gen_fun, df1=10, df2=5) |>
    Spower(replications=30000)
#> 
#> Execution time (H:M:S): 00:00:09
#> Design conditions: 
#> 
#> # A tibble: 1 × 6
#>       n     d   df1   df2 sig.level power
#>   <dbl> <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1   100   0.5    10     5      0.05 NA   
#> 
#> Estimate of power: 0.959
#> 95% Confidence Interval: [0.957, 0.961]

  # evaluate Type I error rate to see if liberal/conservative given
  # assumption violations (should be close to alpha/sig.level)
  p_t.test(n=100, d=0, gen_fun=my.gen_fun, df1=10, df2=5) |>
    Spower(replications=30000)
#> 
#> Execution time (H:M:S): 00:00:09
#> Design conditions: 
#> 
#> # A tibble: 1 × 6
#>       n     d   df1   df2 sig.level power
#>   <dbl> <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1   100     0    10     5      0.05 NA   
#> 
#> Estimate of power: 0.051
#> 95% Confidence Interval: [0.049, 0.054]

# }
```

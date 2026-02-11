# p-value from comparing two or more correlations simulation

Function utilizes [`cocor`](https://rdrr.io/pkg/cocor/man/cocor.html) to
perform correlation comparison for independent, overlapping, and
non-overlapping correlation designs. Type type of correlation design is
inferred based on which correlations are specified.

## Usage

``` r
p_2r(
  n,
  r.ab,
  r.ab2 = NULL,
  r.ac,
  r.bc,
  r.ad,
  r.bd,
  r.cd,
  n2_n1 = 1,
  two.tailed = TRUE,
  test = NULL,
  gen_fun = gen_2r,
  return_analysis = FALSE,
  ...
)

gen_2r(n, R, ...)
```

## Arguments

- n:

  sample size

- r.ab:

  correlation between variable A and B (for independent groups, this is
  for sample 1)

- r.ab2:

  (for independent group test only) correlation between variable A and B
  in sample 2

- r.ac:

  (for overlap/non-overlap) correlation between A and C. This is the
  correlation used in the overlapping hypothesis test, comparing this
  input to `r.ab`

- r.bc:

  (for overlap/non-overlap only) correlation between B and C.

- r.ad:

  (for non-overlap only) correlation between A and D

- r.bd:

  (for non-overlap only) correlation between B and D

- r.cd:

  (for non-overlap only) correlation between C and D. This is the
  correlation used in the non-overlapping hypothesis test, comparing
  this input to `r.ab`

- n2_n1:

  sample size ratio. Only used for independent group test

- two.tailed:

  logical; use two-tailed test?

- test:

  hypothesis testing method to use. Defaults to `'fisher1925'` for the
  independent groups test and `'steiger1980'` for overlap/non-overlap
  tests

- gen_fun:

  function used to generate the required discrete data. Object returned
  must be a `matrix` with `n` rows. Default uses `gen_2r`. User defined
  version of this function must include the argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

- R:

  a correlation matrix constructed from the inputs to `p_2r`

## Value

a single p-value

## Details

For independent group tests, only `r.ab` and `r.ab2` need to be
specified, where the null hypothesis pertains to \\H_0:
r\_{ab}=r\_{ab2}\\.

For overlapping correlation tests, `r.ab`, `r.ac`, and `r.bc` need to be
specified, where the null hypothesis pertains to \\H_0:
r\_{ab}=r\_{ac}\\.

For non-overlapping correlation tests, all correlations expect for
`r.ab2` must be specified, where the null hypothesis pertains to \\H_0:
r\_{ab}=r\_{cd}\\.

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# independent (same x-y pairing across groups)
p_2r(100, r.ab=.5, r.ab2=.6)
#> [1] 0.07403686

# return cocor object for further analysis
p_2r(100, r.ab=.5, r.ab2=.6, return_analysis = TRUE)
#> 
#>   Results of a comparison of two correlations based on independent groups
#> 
#> Comparison between r1.jk (y, x) = 0.5414 and r2.hm (y, x) = 0.5891
#> Difference: r1.jk - r2.hm = -0.0477
#> Data: sample1: j = y, k = x; sample2: h = y, m = x
#> Group sizes: n1 = 100, n2 = 100
#> Null hypothesis: r1.jk is equal to r2.hm
#> Alternative hypothesis: r1.jk is not equal to r2.hm (two-sided)
#> Alpha: 0.05
#> 
#> fisher1925: Fisher's z (1925)
#>   z = -0.4890, p-value = 0.6249
#>   Null hypothesis retained
#> 

# \donttest{

   # estimate empirical power
   p_2r(n=100, r.ab=.5, r.ab2=.6) |> Spower()
#> 
#> Execution time (H:M:S): 00:00:19
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n  r.ab r.ab2 sig.level power
#>   <dbl> <dbl> <dbl>     <dbl> <lgl>
#> 1   100   0.5   0.6      0.05 NA   
#> 
#> Estimate of power: 0.170
#> 95% Confidence Interval: [0.163, 0.177]

   # estimate n required to reach 80% power
   p_2r(n=NA, r.ab=.5, r.ab2=.6) |>
        Spower(power=.80, interval=c(100, 5000))
#> 
#> Execution time (H:M:S): 00:01:19
#> Design conditions: 
#> 
#> # A tibble: 1 × 5
#>       n  r.ab r.ab2 sig.level power
#>   <dbl> <dbl> <dbl>     <dbl> <dbl>
#> 1    NA   0.5   0.6      0.05   0.8
#> 
#> Estimate of n: 762.4
#> 95% Predicted Confidence Interval: [752.9, 771.4]

# }

# overlap (same y, different xs)
# H0: r.ab = r.bc
p_2r(100, r.ab=.5, r.ac=.3, r.bc=.2)
#> [1] 0.01196608

# nonoverlap (different ys, different xs)
# H0: r.ab = r.cd
p_2r(100, r.ab=.5, r.ac=.3, r.bc=.2, r.ad=.2, r.bd=.4, r.cd=.2)
#> [1] 0.7617057

```

# Evaluate whether a confidence interval is within a tolerable interval

Return `TRUE` if an estimated confidence interval falls within a
tolerable `interval` range. Typically used for equivalence, superiority,
or non-inferiority testing.

## Usage

``` r
is.CI_within(CI, interval)
```

## Arguments

- CI:

  estimated confidence interval (length 2)

- interval:

  tolerable interval range (length 2)

## Value

logical

## See also

[`is.outside_CI`](https://philchalmers.github.io/Spower/reference/is.outside_CI.md),
[`Spower`](https://philchalmers.github.io/Spower/reference/Spower.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
CI <- c(.2, .4)
LU <- c(.1, .3)
is.CI_within(CI, LU)        # not within tolerable interval
#> [1] FALSE
is.CI_within(CI, c(0, .5))  # is within wider interval
#> [1] TRUE

# complement indicates if CI is outside interval
!is.CI_within(CI, LU)
#> [1] TRUE

#####
# for superiority test
is.CI_within(CI, c(.1, Inf))  # CI is within tolerable interval
#> [1] TRUE

# for inferiority test
is.CI_within(CI, c(-Inf, .3))  # CI is not within tolerable interval
#> [1] FALSE
```

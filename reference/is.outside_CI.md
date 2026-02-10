# Evaluate whether parameter is outside a given confidence interval

Returns `TRUE` if parameter reflecting a null hypothesis falls outside a
given confidence interval. This is an alternative approach to writing an
experiment that returns a p-value.

## Usage

``` r
is.outside_CI(P0, CI)
```

## Arguments

- P0:

  parameter to evaluate

- CI:

  confidence interval

## Value

logical

## See also

[`is.CI_within`](https://philchalmers.github.io/Spower/reference/is.CI_within.md),
[`Spower`](https://philchalmers.github.io/Spower/reference/Spower.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
p0 <- .3
CI <- c(.2, .4)
is.outside_CI(p0, CI)
#> [1] FALSE

# complement indicates if p0 is within CI
!is.outside_CI(p0, CI)
#> [1] TRUE

```

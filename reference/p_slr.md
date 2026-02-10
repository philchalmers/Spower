# p-value from simple linear regression model simulation

p-values associated with the simple linear regression model, \\y =
\beta_0 + \beta_1 X + \epsilon\\. Focus is on the slope/intercept
behavior of the model.

## Usage

``` r
p_slr(
  n,
  beta,
  sd_x,
  sd_y,
  beta0 = 0,
  test = "x = 0",
  gen_fun = gen_slr,
  return_analysis = FALSE,
  ...
)

gen_slr(n, beta, sd_x, sd_y, ...)
```

## Arguments

- n:

  sample size

- beta:

  slope parameter

- sd_x:

  standard deviation of IV

- sd_y:

  standard deviation of DV

- beta0:

  null value to test against

- test:

  test to evaluate using
  [`lht`](https://rdrr.io/pkg/car/man/linearHypothesis.html). Default
  evaluates the null hypothesis that the slope is equal to 0

- gen_fun:

  function used to generate the required X-Y data. Object returned must
  be a `data.frame` with the columns `'y'` and `'x'`. Default uses
  `gen_slr`. User defined version of this function must include the
  argument `...`

- return_analysis:

  logical; return the analysis object for further extraction and
  customization?

- ...:

  additional arguments to be passed to `gen_fun`. Not used unless a
  customized `gen_fun` is defined

## Value

a single p-value

## See also

[`p_glm`](https://philchalmers.github.io/Spower/reference/p_glm.md),
[`p_lm.R2`](https://philchalmers.github.io/Spower/reference/p_lm.R2.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
p_slr(n=100, beta = -0.0667, sd_x = 7.5, sd_y = 4)
#> [1] 0.5339516

if (FALSE) { # \dontrun{
p_slr(n=100, beta = -0.0667, sd_x = 7.5, sd_y = 4) |> Spower()
} # }
```

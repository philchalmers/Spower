
[![](http://www.r-pkg.org/badges/version/Spower)](https://www.r-pkg.org:443/pkg/Spower)
[![](http://cranlogs.r-pkg.org/badges/grand-total/Spower)](https://CRAN.R-project.org/package=Spower)

# *Spower*: Power Analyses using Monte Carlo Simulations <img src="inst/sticker/S.png" height="139" align="right"/>

*Spower* provides a general purpose simulation-based power analysis API
for routine and customized simulation experimental designs. The package
focuses exclusively on Monte Carlo simulation variants of (expected)
prospective power analyses, a priori power analyses, criterion power
analyses, compromise power analyses, and sensitivity analyses. The
default simulation experiment functions found within the package provide
stochastic variants of the power analyses subroutines found in the
*GPower* 3 software (Faul, Erdfelder, Buchner, and Lang, 2009) along
with various other power analysis examples (e.g., mediation analyses).
Supporting functions are also included, such as for building empirical
power curve estimates, which utilize a similar API structure.

## Installation

You can install `Spower` from CRAN:

``` r
install.packages("Spower")
```

To install the development version of the `Spower` package, you need to
install the `remotes` package then the `Spower` package.

``` r
install.packages("remotes")
remotes::install_github("philchalmers/Spower")
```

## Quick Start

*Spower* requires only two components: an available function used to
generate exactly one simulation experiment that returns one or more
*p*-values given the null hypothesis of interest, and the use of either
`Spower()` or `SpowerCurve()` to perform the desired
prospective/post-hoc, a priori, sensitivity, compromise, or criterion
power analysis.

For example, the built-in `p_t.test()` function performs *t*-tests using
various inputs, where below a sample size of $N=200$ is supplied
(`n = 100` per group) and a Cohen’s $d$ of .5 (medium effect). This
returns a single $p$-value given the null hypotheses tested of no mean
difference, which in this single case returns a ‘surprising’ result
given this null position tested.

``` r
library(Spower)
p_t.test(n=100, d=0.5)
## [1] 0.001231514
```

To evaluate the prospective power of this test requires passing the
simulation function to `Spower()`.

``` r
set.seed(42)
p_t.test(n=100, d=0.5) |> Spower()

##
## Execution time (H:M:S): 00:00:02
## Design conditions: 
##
## # A tibble: 1 × 4
##       n     d sig.level power
##   <dbl> <dbl>     <dbl> <lgl>
## 1   100   0.5      0.05 NA   
## 
## Estimate of power: 0.943
## 95% Confidence Interval: [0.938, 0.947]
```

To evaluate the prospective power of this test requires passing the
simulation function to `Spower()`. For a priori and sensitive analyses,
the respective input to the simulation function must be set to `NA`,
while within `Spower()` the target power rate must be included along
with a suitable search `interval`.

``` r
set.seed(01123581321)

# estimate the require n value to achieve a power of 1 - beta = .95 
p_t.test(n=NA, d=0.5) |> Spower(power=.95, interval=c(50, 300))

## Iter: 56; Median = 101; E(f(x)) = 0.00; Total.reps = 11000; k.tol = 2; Pred = 103.4
##
## Execution time (H:M:S): 00:00:05
## Design conditions: 
## 
## # A tibble: 1 × 4
##       n     d sig.level power
##   <dbl> <dbl>     <dbl> <dbl>
## 1    NA   0.5      0.05  0.95
## 
## Estimate of n: 103.4
## 95% Predicted Confidence Interval: [101.5, 105.4]
```

To generate suitable power-curves for any given simulation or power
analysis criteria the simulation experiment can be passed to
`SpowerCurve()`. This function contains a similar specification
structure to `Spower()`, however differs in that the arguments to vary
are explicitly passed as named vectors to `SpowerCurve()`. Below is an
example varying sample size (`n`), while the next example varies both
`n` and the effect size `d`.

``` r
set.seed(8675309) # Jenny Jenny
p_t.test(d=0.2) |> SpowerCurve(n=seq(50, 550, by=100))
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

``` r
p_t.test() |> SpowerCurve(n=seq(50, 350, by=50), d=c(.2, .3, .4, .5))
```

![](man/figures/README-unnamed-chunk-6-1.png)<!-- -->

## Vignettes

The package currently contains a vignette demonstrating several of the
examples from the *GPower 3.1* manual, providing simulation-based
replications of the same analyses, as well other vignettes showing more
advanced applications of the package (ROPEs, Bayesian power anlayses,
precision criterion, etc).

## Getting Help or Reporting an Issue

To report bugs/issues/feature requests, please file an
[issue](https://github.com/philchalmers/Spower/issues/).

## How to Contribute

If you have a simulation experiment you’d like to contribute in the form
of

``` r
# returns a p-value, P(D|H0)
p_yourSimulation()

# returns a posterior probability, P(H1|D)
pp_yourSimulation()

# returns a logical (more complex experiments, perhaps involing ROPEs)
l_yourSimulation()
```

then feel free to document this function using the `roxygen2` style
syntax and open a pull request. Otherwise, contributions can be made to
the online Wiki.

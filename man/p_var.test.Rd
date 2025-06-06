% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p_var.test.R
\name{p_var.test}
\alias{p_var.test}
\alias{gen_var.test}
\title{p-value from variance test simulation}
\usage{
p_var.test(
  n,
  vars,
  n.ratios = rep(1, length(vars)),
  sigma2 = 1,
  two.tailed = TRUE,
  test = "Levene",
  correct = TRUE,
  gen_fun = gen_var.test,
  ...
)

gen_var.test(n, vars, n.ratios = rep(1, length(vars)), ...)
}
\arguments{
\item{n}{sample size per group, assumed equal across groups}

\item{vars}{a vector of variances to use for each group; length of 1 for
one-sample tests}

\item{n.ratios}{allocation ratios reflecting the sample size ratios.
Default of 1 sets the groups to be the same size (n * n.ratio)}

\item{sigma2}{population variance to test against in one-sample test}

\item{two.tailed}{logical; should a two-tailed or one-tailed test be used?}

\item{test}{type of test to use in multi-sample applications.
Can be either \code{'Levene'} (default), \code{'Bartlett'}, or \code{'Fligner'}}

\item{correct}{logical; use correction when \code{test = 'Bartlett'}?}

\item{gen_fun}{function used to generate the required discrete data.
Object returned must be a \code{matrix} with k rows and k columns
of counts. Default uses \code{\link{gen_var.test}}.
User defined version of this function must include the argument \code{...}}

\item{...}{additional arguments to be passed to \code{gen_fun}. Not used
unless a customized \code{gen_fun} is defined}
}
\value{
a single p-value
}
\description{
Generates one or or more sets of continuous data group-level data
to perform a variance test, and return a p-value. When two-samples
are investigated the \code{\link{var.test}} function will be used,
otherwise functions from the \code{EnvStats} package will be used.
}
\examples{

# one sample
p_var.test(100, vars=10, sigma2=9)

# three sample
p_var.test(100, vars=c(10, 9, 11))
p_var.test(100, vars=c(10, 9, 11), test = 'Fligner')
p_var.test(100, vars=c(10, 9, 11), test = 'Bartlett')

\donttest{
  # power to detect three-group variance differences
  p_var.test(n=100, vars=c(10,9,11)) |> Spower()

  # sample size per group to achieve 80\% power
  p_var.test(n=NA, vars=c(10,9,11)) |>
         Spower(power=.80, interval=c(100, 1000))
}

}
\seealso{
\code{\link{gen_var.test}}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}

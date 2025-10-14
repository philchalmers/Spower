#' p-value from Kolmogorov-Smirnov one- or two-sample simulation
#'
#' Generates one or two sets of continuous data group-level data
#' and returns a p-value under the null that the groups were drawn from
#' the same distribution (two sample) or from a theoretically known distribution
#' (one sample).
#'
#' @param n sample size per group, assumed equal across groups
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param p1 a function indicating how the data were generated for group 1
#' @param p2 (optional) a function indicating how the data were generated for group 2.
#'   If omitted a one-sample test will be evaluated provided that \code{parent} is also
#'   specified
#' @param n2_n1 sample size ratio. Default uses equal sample sizes
#' @param parent the cumulative distribution function to use
#'   (e.g., \code{\link{pnorm}}). Specifying this input will construct a
#'   one-sample test setup
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @param ... additional arguments to be passed to the
#'   \code{parent} distribution function from \code{\link{ks.test}}, as well as any other
#'   relevant parameter to \code{ks.test} (e.g., \code{exact = TRUE})
#' @seealso \code{\link{gen_t.test}}
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @examples
#'
#' # two-sample test from two Gaussian distributions with different locations
#' p1 <- function(n) rnorm(n)
#' p2 <- function(n) rnorm(n, mean=-.5)
#' p_ks.test(n=100, p1, p2)
#'
#' # return analysis model
#' p_ks.test(n=100, p1, p2, return_analysis=TRUE)
#'
#' # one-sample data from chi-squared distribution tested
#' #   against a standard normal distribution
#' pc <- function(n, df=15) (rchisq(n, df=df) - df) / sqrt(2*df)
#' p_ks.test(n=100, p1=pc, parent=pnorm, mean=0, sd=1)
#'
#' \donttest{
#'   # empirical power estimates
#'   p_ks.test(n=100, p1, p2) |> Spower()
#'   p_ks.test(n=100, p1=pc, parent=pnorm, mean=0, sd=1) |> Spower()
#'
#' }
#'
#' @export
p_ks.test <- function(n, p1, p2, n2_n1 = 1, two.tailed = TRUE,
					  parent = NULL, return_analysis = FALSE, ...) {
	n2 <- n2_n1 * n
	dat1 <- p1(n=n)
	if(is.null(parent)){
		dat2 <- p2(n=n2)
	} else dat2 <- parent
	ret <- ks.test(dat1, dat2, ...)
	if(return_analysis) return(ret)
	p <- ret$p.value
	p <- ifelse(two.tailed, p, p/2)
	p
}


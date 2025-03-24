#' p-value from Wilcox test simulation
#'
#' Simulates data given one or two parent distributions and
#' returns a p-value. Can also be used for power analyses related
#' to sign tests.
#'
#' @param n sample size per group
#' @param d effect size passed to \code{parent} functions
#' @param n2_n1 sample size ratio
#' @param type type of analysis to use (two-sample, one-sample, or paired)
#' @param mu parameter used to form the null hypothesis
#' @param exact a logical indicating whether an exact p-value should be computed
#' @param correct a logical indicating whether to apply continuity correction
#'   in the normal approximation for the p-value
#' @param parent1 data generation function for first group. Ideally
#'   should have SDs = 1 so that \code{d} reflects a standardized
#'   difference
#' @param parent2 same as \code{parent1}, but for the second group
#' @param two.tailed logical; use two-tailed test?
#' @export
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @examples
#'
#' # with normal distributions defaults d is standardized
#' p_wilcox.test(100, .5)
#' p_wilcox.test(100, .5, type = 'paired')
#' p_wilcox.test(100, .5, type = 'one.sample')
#'
#' # using chi-squared distributions (standardizing to 0-1)
#' p_wilcox.test(100, .5, type = 'one.sample',
#'    parent1 = function(n, d) rchisq(n, df=10) - 10 + d)
#' p_wilcox.test(100, .5,
#'    parent1 = function(n, d) (rchisq(n, df=10) - 10)/sqrt(20) + d,
#'    parent2 = function(n, d) (rchisq(n, df=10) - 10)/sqrt(20))
#'
p_wilcox.test <- function(n, d, n2_n1 = 1, mu=0,
						  type = c('two.sample', 'one.sample', 'paired'),
						  exact = NULL, correct = TRUE, two.tailed = TRUE,
						  parent1 = function(n, d) rnorm(n, d, 1),
						  parent2 = function(n, d) rnorm(n, 0, 1)){
	type <- match.arg(type)
	if(type == 'paired') n <- n * 2
	dat1 <- parent1(n, d)
	ret <- if(type == 'one.sample'){
		wilcox.test(dat1, mu=mu, correct=correct, exact=exact)
	} else {
		dat2 <- parent2(n*n2_n1, d)
		paired <- type == 'paired'
		wilcox.test(dat1, dat2, paired=paired,
					mu=mu, correct=correct, exact=exact)
	}
	p <- ret$p.value
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' p-value from Wilcoxon (signed rank) test simulation
#'
#' Simulates data given one (Wilcoxon) or two (Mann-Whitney)
#' parent distributions and returns a p-value. Can also be used for power analyses related
#' to sign tests.
#'
#' @param n sample size per group. For paired samples this corresponds to the
#'   number of pairings (hence, half of the data points observed)
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
#'   For one-tailed tests relevant \code{alternative} arguments are set to \code{'greater'}
#' @param alternative character vector passed to \code{alternative} analysis function
#'   (most commonly, \code{'two.tailed', 'greater', 'less'}). When used, ignores the value
#'   of \code{two.tailed}
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @export
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @examples
#'
#' # with normal distributions defaults d is standardized
#' p_wilcox.test(100, .5)
#' p_wilcox.test(100, .5, type = 'paired')  # n = number of pairs
#' p_wilcox.test(100, .5, type = 'one.sample')
#'
#' # return analysis object
#' p_wilcox.test(100, .5, return_analysis = TRUE)
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
						  alternative = NULL,
						  parent1 = function(n, d) rnorm(n, d, 1),
						  parent2 = function(n, d) rnorm(n, 0, 1),
						  return_analysis = FALSE){
	type <- match.arg(type)
	dat1 <- parent1(n, d)
	if(is.null(alternative))
		alternative <- ifelse(two.tailed, "two.sided", "greater")
	ret <- if(type == 'one.sample'){
		wilcox.test(dat1, mu=mu, correct=correct, exact=exact,
					alternative=alternative)
	} else {
		paired <- type == 'paired'
		if(paired) stopifnot(n2_n1 == 1)
		dat2 <- parent2(n*n2_n1, d)
		wilcox.test(dat1, dat2, paired=paired, alternative=alternative,
					mu=mu, correct=correct, exact=exact)
	}
	if(return_analysis) return(ret)
	p <- ret$p.value
	p
}

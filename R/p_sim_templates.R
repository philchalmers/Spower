#' Independent samples t-test simulation and p-value
#'
#' Function generates two sets of individuals with equal sample sizes
#' and Cohen's effect size 'd'. The data and associated t-test
#' assume that the conditional observations are normally distributed and have
#' have equal variance.
#'
#' @param n sample size per group, assumed equal across groups
#' @param d Cohen's standardized effect size
#' @return a single p-value
#' @examples
#'
#' # sample size of 50 per group, medium effect size
#' p_t.test(n=50, d=0.5)
#'
#' @export
p_t.test <- function(n, d) {
	group1 <- rnorm(n)
	group2 <- rnorm(n, mean=d)
	dat <- data.frame(group = gl(2, n, labels=c('G1', 'G2')),
					  DV = c(group1, group2))
	p <- t.test(DV ~ group, dat, var.equal=TRUE)$p.value
	p
}

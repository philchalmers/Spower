#' Independent/paired samples t-test simulation and p-value
#'
#' Function generates two sets of individuals with equal sample sizes
#' and Cohen's effect size 'd'. The data and associated t-test
#' assume that the conditional observations are normally distributed and have
#' have equal variance.
#'
#' @param n sample size per group, assumed equal across groups
#' @param d Cohen's standardized effect size
#' @param mu population mean to test against
#' @param paired logical; should the analysis results assume that the groups
#'   are dependent?
#' @param n2_n1 allocation ratio reflecting the same size ratio.
#'   Default of 1 sets the groups to be the same size. Only applicable
#'   when \code{paired = FALSE}
#' @return a single p-value
#' @examples
#'
#' # sample size of 50 per group, medium effect size
#' p_t.test(n=50, d=0.5)
#' p_t.test(n=50, d=0.5, paired=TRUE)
#'
#' # second group 2x as large as the first group
#' p_t.test(n=50, d=0.5, n2_n1 = 2)
#'
#' @export
p_t.test <- function(n, d, mu = 0, paired = FALSE,
					 n2_n1 = 1, two.tailed = TRUE) {
	group1 <- rnorm(n)
	group2 <- rnorm(n * n2_n1, mean=d)
	dat <- data.frame(group = factor(rep(c('G1', 'G2'), times=c(n, n*n2_n1))),
					  DV = c(group1, group2))
	alternative <- ifelse(two.tailed, 'two.sided',
						  ifelse(d < 0, 'greater', 'less'))
	p <- if(paired){
		if(n2_n1 != 1) stop('n2_n1 must equal 1 for paired t-tests')
		t.test(group1, group2, mu=mu, paired=TRUE,
			   alternative=alternative)$p.value
	} else t.test(DV ~ group, dat, var.equal=TRUE,
				  alternative=alternative, mu=mu)$p.value
	p
}

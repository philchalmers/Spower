#' Independent/paired samples t-test simulation and p-value
#'
#' Function generates two sets of individuals with equal sample sizes
#' and Cohen's effect size 'd'. The data and associated t-test
#' assume that the conditional observations are normally distributed and have
#' have equal variance.
#'
#' @param n sample size per group, assumed equal across groups
#' @param d Cohen's standardized effect size \code{d}
#' @param mu population mean to test against
#' @param paired logical; should the analysis results assume that the groups
#'   are dependent?
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param n2_n1 allocation ratio reflecting the same size ratio.
#'   Default of 1 sets the groups to be the same size. Only applicable
#'   when \code{paired = FALSE}
#' @param raw_info (optional) list of mean and variance inputs for each group.
#'   If supplied the \code{d} input will be ignored
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
					 n2_n1 = 1, two.tailed = TRUE,
					 raw_info = list(mu1=NA, mu2=NA,
					 				sigma1=NA, sigma2=NA)) {
	if(!is.na(raw_info$mu1)){
		group1 <- with(raw_info, rnorm(n, mean=mu1, sd=sigma1))
		group2 <- with(raw_info, rnorm(n * n2_n1, mean=mu2, sd=sigma2))
	} else {
		group1 <- rnorm(n)
		group2 <- rnorm(n * n2_n1, mean=d)
	}
	dat <- data.frame(group = factor(rep(c('G1', 'G2'), times=c(n, n*n2_n1))),
					  DV = c(group1, group2))
	p <- if(paired){
		if(n2_n1 != 1) stop('n2_n1 must equal 1 for paired t-tests')
		t.test(group1, group2, mu=mu, paired=TRUE)$p.value
	} else t.test(DV ~ group, dat, var.equal=TRUE, mu=mu)$p.value
	p <- ifelse(two.tailed, p, p*2)
	p
}

#' Correlation simulation and p-value
#'
#' Function generates correlated X-Y data and returns a p-value to assess the null
#' of no correlation in the population. The X-Y data are generated
#' assuming a multivariate normal distribution.
#'
#' @param n sample size
#' @param r correlation
#' @param method method to use to compute the correlation
#'   (see \code{\link{cor.test}})
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @return a single p-value
#' @examples
#'
#' # 50 observations, .5 correlation
#' p_r(50, r=.5)
#' p_r(50, r=.5, method = 'spearman')
#'
#' @export
p_r <- function(n, r, rho = 0, method = 'pearson', two.tailed = TRUE) {
	dat <- SimDesign::rmvnorm(n, sigma = matrix(c(1,r,r,1), 2, 2))
	colnames(dat) <- c('x', 'y')
	test <- cor.test(~ x + y, dat, method=method)
	p <- ifelse(two.tailed, test$p.value, test$p.value*2)
	p
}

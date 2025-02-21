#' p-value from one-way ANOVA simulation
#'
#' Generates continuous multi-sample data to be analysed by
#' a one-way ANOVA, and return a p-value.
#' Uses the function \code{\link{oneway.test}} to perform the analyses.
#' The data and associated
#' test assume that the conditional observations are normally distributed and have
#' have equal variance by default, however these may be modified.
#'
#' @param n sample size per group
#' @param k number of groups
#' @param f Cohen's f effect size
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param var.equal logical; use the pooled SE estimate instead of the Welch
#'   correction for unequal variances?
#' @param means (optional) vector of means. When specified the input \code{f} is ignored
#' @param sds (optional) vector of SDs. When specified the input \code{f} is ignored
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{matrix} with k rows and k columns
#'   of counts. Default uses \code{\link{gen_anova.test}}.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#' @return a single p-value
#' @seealso \code{\link{gen_anova.test}}
#' @examples
#'
#' # n=50 in 3 groups, "medium" effect size
#' p_anova.test(50, k=3, f=.25)
#'
#' # explicit means/sds
#' p_anova.test(50, 3, means=c(0,0,1), sds=c(1,2,1))
#'
#' if(FALSE){
#'   # compare simulated results to pwr package
#'   pwr::pwr.anova.test(f=0.28, k=4, n=20)
#'   p_anova.test(n=20, k=4, f=.28) |> Spower()
#' }
#'
#' @export
p_anova.test <- function(n, k, f, n.ratios = rep(1, k),
						 two.tailed = TRUE, var.equal = TRUE,
						 means=NULL, sds=NULL, gen_fun=gen_anova.test, ...) {
	df <- gen_fun(n=n, k=k, f=f, n.ratios=n.ratios, means=means, sds=sds, ...)
	p <- oneway.test(DV ~ group, data=df, var.equal=var.equal, subset=NULL,
					 na.action = NULL)$p.value
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' @rdname p_anova.test
#' @export
gen_anova.test <- function(n, k, f, n.ratios = rep(1, k),
						   means=NULL, sds=NULL, ...){
	stopifnot(length(n) == 1)
	stopifnot(length(n.ratios) == k)
	group <- rep(factor(1:k), times = n*n.ratios)
	n.each <- n*n.ratios
	stopifnot(all.equal(n.each, as.integer(n.each)))
	if(!is.null(means)){
		dv <- sapply(1:k, \(i)
					 rnorm(n*n.ratios[i], mean=means[i], sd=sds[i]))
		df <- data.frame(group=group, DV = as.numeric(dv))
	} else {
		# use f
		f2 <- f^2
		N <- sum(n.each)
		eta2 <- 1 / (1/f2 + 1)
		# pick a nice SSG from mu=0 deviation form, and solve for SSE
		gmeans <- seq(-2, 2, length.out=k)
		SSG <- sum(n.each * gmeans^2)
		SSE <- SSG*(1/eta2 - 1)
		sde <- sqrt(SSE / N)
		dv <- rep(gmeans, times=n.each) + rnorm(N, sd=sde)
		df <- data.frame(DV=dv, group=group)
	}
	df
}

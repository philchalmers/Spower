#' p-value from Kruskal-Wallis Rank Sum Test simulation
#'
#' Simulates data given two or more parent distributions and
#' returns a p-value using \code{\link{kruskal.test}}. Default generates data
#' from Gaussian distributions, however this can be modified.
#'
#' @param n sample size per group
#' @param k number of groups
#' @param means vector of means to control location parameters
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param gen_fun function used to generate the required data.
#'   Object returned must be a \code{list} of length \code{k}, where each element
#'   contains the sample data in each group. Default uses \code{\link{gen_kruskal.test}}.
#'   User defined version of this function must include the argument \code{...}
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @param ... additional arguments to pass to \code{gen_fun}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @return a single p-value
#' @export
#' @examples
#'
#' # three group test where data generate from Gaussian distributions
#' p_kruskal.test(n=30, k=3, means=c(0, .5, .6))
#'
#' # return analysis model
#' p_kruskal.test(n=30, k=3, means=c(0, .5, .6), return_analysis=TRUE)
#'
#' # generate data from chi-squared distributions with different variances
#' gen_chisq <- function(n, k, n.ratios, means, dfs, ...){
#'   dat <- vector('list', k)
#'   ns <- n * n.ratios
#'   for(g in 1:k)
#'  	dat[[g]] <- rchisq(ns[g], df=dfs[g]) - dfs[g] + means[g]
#'   dat
#' }
#'
#' p_kruskal.test(n=30, k=3, means=c(0, 1, 2),
#'    gen_fun=gen_chisq, dfs=c(10, 15, 20))
#'
#' \donttest{
#'   # empirical power estimate
#'   p_kruskal.test(n=30, k=3, means=c(0, .5, .6)) |> Spower()
#'   p_kruskal.test(n=30, k=3, means=c(0, 1, 2), gen_fun=gen_chisq,
#'          dfs = c(10, 15, 20)) |> Spower()
#'
#' }
#'
p_kruskal.test <- function(n, k, means, n.ratios = rep(1, k),
							gen_fun=gen_kruskal.test, return_analysis = FALSE, ...){
	dat <- gen_fun(n, k=k, n.ratios=n.ratios, means=means, ...)
	ret <- kruskal.test(dat)
	if(return_analysis) return(ret)
	p <- ret$p.value
	p
}

#' @rdname p_kruskal.test
#' @export
gen_kruskal.test <- function(n, k, n.ratios, means, ...){
	dat <- vector('list', k)
	ns <- n * n.ratios
	for(g in 1:k)
		dat[[g]] <- rnorm(ns[g], mean=means[g])
	dat
}

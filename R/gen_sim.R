#' Generate data for paired, one-sample, and two-sample t.test procedures
#'
#' Generates one or two sets of continuous data group-level data
#' according to Cohen's effect size 'd'. The data are generated such that
#' the conditional observations are normally distributed and have
#' have equal variance, however some of these properties may be modified.
#'
#' @param n sample size per group, assumed equal across groups
#' @param d Cohen's standardized effect size \code{d}
#' @param r (optional) instead of specifying \code{d} specify
#'   a point-biserial correlation. Internally this is transformed
#'   into a suitable \code{d} value for the power computations
#' @param type type of data structure to generate; can be \code{'two.sample'},
#'   \code{'one.sample'}, or \code{'paired'}
#' @param n2_n1 allocation ratio reflecting the same size ratio.
#'   Default of 1 sets the groups to be the same size. Only applicable
#'   when \code{type = 'two.sample'}
#' @param raw_info (optional) list of mean and SD inputs for each group,
#'   each specified as a vector. When specified the input \code{d} is ignored
#' @param ... additional arguments (not used)
#'
#' @return a \code{data.frame} with the columns \code{'DV'} and \code{'group'}
#'
#' @seealso \code{\link{p_t.test}}
#' @export
#' @examples
#'
#' # two-sample data with 100 observations each
#' df <- gen_t.test(n=100, d=.5)
#' head(df)
#' with(df, tapply(DV, group, mean))
#'
gen_t.test <- function(n, d, n2_n1 = 1, r = NULL,
					   type = c('two.sample', 'one.sample', 'paired'),
					   raw_info = list(means=NA, sds=NA), ...){
	type <- match.arg(type)
	if(!is.null(r)){
		type <- 'two.sample'
		stopifnot(n2_n1 == 1)
	}
	if(!is.null(r))
		d <- r2d(r, n0=n, n1=n*n2_n1)
	if(type == 'paired') n <- n * 2
	n.each <- n * n2_n1
	stopifnot(all.equal(n.each, as.integer(n.each)))
	if(type == 'one.sample'){
		DV <- if(!all(is.na(raw_info$means)))
			with(raw_info, rnorm(n, mean=means, sd=sds)) else rnorm(n, mean=d)
		dat <- data.frame(DV=DV)
	} else {
		if(!all(is.na(raw_info$means))){
			if(!missing(d)) stop('d argument cannot be used with raw_info')
			group1 <- with(raw_info, rnorm(n, mean=means[1], sd=sds[1]))
			group2 <- with(raw_info, rnorm(n * n2_n1, mean=means[2], sd=sds[2]))
		} else {
			group1 <- rnorm(n)
			group2 <- rnorm(n * n2_n1, mean=d)
		}
		dat <- data.frame(group = factor(rep(c('G1', 'G2'), times=c(n, n*n2_n1))),
						  DV = c(group1, group2))
	}
	dat
}

#' Generate bivariate normal data
#'
#' Generates correlated X-Y data from a bivariate normal distribution.
#'
#' @param n sample size
#' @param r correlation
#' @param ... additional arguments (not used)
#' @return a \code{matrix} with two-columns and \code{n} rows
#' @export
#' @seealso \code{\link{p_r}} and \code{\link{p_r.cat}}
#' @examples
#'
#' dat <- gen_r(1000, r=.9)
#' plot(dat)
#'
gen_r <- function(n, r, ...){
	dat <- SimDesign::rmvnorm(n, sigma = matrix(c(1,r,r,1), 2, 2))
	dat
}

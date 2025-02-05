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
#' @param means (optional) vector of means for each group.
#'   When specified the input \code{d} is ignored
#' @param sds (optional) vector of SDs for each group.
#'   When specified the input \code{d} is ignored
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
					   means=NULL, sds=NULL, ...){
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
		DV <- if(!is.null(means))
			rnorm(n, mean=means, sd=sds) else rnorm(n, mean=d)
		dat <- data.frame(DV=DV)
	} else {
		if(!is.null(means)){
			if(!missing(d)) stop('d argument cannot be used with raw_info')
			group1 <- rnorm(n, mean=means[1], sd=sds[1])
			group2 <- rnorm(n * n2_n1, mean=means[2], sd=sds[2])
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

#' Generate multivariate normal data
#'
#' Generates correlated X-Y data from a bivariate normal distribution.
#'
#' @param n sample size
#' @param mean vector of means. Defaults to a vector of 0's
#' @param sigma covariance matrix, though if all diagonal values equal 1
#'   then understood as a correlation matrix
#' @param ... additional arguments (not used)
#' @return a \code{matrix} with two-columns and \code{n} rows
#' @export
#' @seealso \code{\link{p_r}} and \code{\link{p_r.cat}}
#' @examples
#'
#' sigma <- matrix(c(1,.2,.3,.2,1,.5,.3,.5,1), 3, 3)
#' dat <- gen_mvnorm(1000, sigma=sigma)
#' pairs(dat)
#'
gen_mvnorm <- function(n, mean = numeric(nrow(sigma)), sigma, ...){
	dat <- SimDesign::rmvnorm(n, mean=mean, sigma=sigma)
	dat
}

#' Generate binomial data for proportion tests
#'
#' Generates single and multi-sample binomial data for proportion tests.
#'
#' @param n sample size per group
#' @param prop sample probability/proportions of success.
#'   If a vector with two-values or more elements are supplied then
#'   a multi-samples test will be used. Matrices are also supported
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param h Cohen's h effect size; only supported for one-sample analysis
#' @param pi probability of success to test against (default is .5). Ignored
#'   for two-sample tests
#' @param ... additional arguments (not used)
#' @return an integer vector for one-sample data generation or
#'   a 2xk matrix of counts for multi-sample problems
#' @seealso \code{\link{p_prop.test}}
#' @examples
#'
#' # one sample, 50 observations
#' gen_prop.test(50, prop=.65)
#'
#' # specified using h and pi
#' h <- pwr::ES.h(.65, .4)
#' gen_prop.test(50, h=h, pi=.4)
#' gen_prop.test(50, h=-h, pi=.65)
#'
#' # two-sample test
#' gen_prop.test(50, prop=c(.5, .65))
#'
#' # two-sample test, unequal ns
#' gen_prop.test(50, prop=c(.5, .65), n.ratios = c(1,2))
#'
#' # three-sample test, group2 twice as large as others
#' gen_prop.test(50, prop=c(.5, .65, .7), n.ratios=c(1,2,1))
#'
#' # data for Fisher exact test
#' gen_prop.test(50, prop=matrix(c(.5, .65, .7, .5), 2, 2))
#'
#' @export
gen_prop.test <- function(n, h, prop, pi = .5, n.ratios = rep(1, length(prop)), ...) {
	if(!missing(h)){
		root.h <- function(p1, p2, h)
			(2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))) - h
		n.ratios <- 1
		int <- if(h > 0) c(pi, 1) else c(0, pi)
		root <- uniroot(root.h, interval=int,
						p2=ifelse(length(n.ratios) == 2, prop, pi), h=h)$root
		prop <- root
		# pwr::ES.h(prop, pi) == h
	}
	stopifnot(length(n) == 1)
	n.each <- n * n.ratios
	stopifnot(all.equal(n.each, as.integer(n.each)))
	dat <- if(length(prop) > 1){
		sapply(1:length(prop), \(i){
			vals <- rbinom(n * n.ratios[i], 1, prob=prop[i])
			c(sum(vals), length(vals))
		})
	} else {
		matrix(c(sum(rbinom(n, 1, prob = prop)), n), nrow=2)
	}
	dat
}

#' Generate multinomial data for McNemar test
#'
#' Generates two-dimensional sample data for McNemar test
#' from a multinomial distribution.
#'
#' @param n total sample size
#' @param prop two-dimensional matrix of proportions/probabilities
#' @param ... additional arguments (not used)
#' @return a matrix of counts
#' @seealso \code{\link{p_mcnemar.test}}
#' @examples
#'
#' # from ?mcnemar.test
#' Performance <- matrix(c(794, 86, 150, 570),
#' 		   nrow = 2,
#' 		   dimnames = list("1st Survey" = c("Approve", "Disapprove"),
#' 		               "2nd Survey" = c("Approve", "Disapprove")))
#' (prop <- prop.table(Performance))
#'
#' # table of counts given prop input
#' gen_mcnemar.test(n=sum(Performance), prop=prop)
#'
#'
#' @export
gen_mcnemar.test <- function(n, prop, ...) {
	draws <- rmultinom(1, n, prob = as.numeric(prop))
	dat <- matrix(draws, nrow(prop), ncol(prop))
	dat
}

#' Generate multinomial data for chi-squared test
#'
#' Generates multinomial data suitable for analysis with
#' \code{\link{chisq.test}}.
#'
#' @param n sample size per group
#' @param P specific power configuration, specified as a numeric vector or matrix
#' @param ... additional arguments (not used)
#' @export
#' @return a numeric vector or matrix, depending on the supplied \code{P} class
#' @seealso \code{\link{p_chisq.test}}
#' @examples
#'
#' # vector of explicit probabilities (goodness of fit test)
#' gen_chisq.test(100, P0 = c(.25, .25, .25, .25),
#'                      P = c(.6, .2, .1, .1))
#'
#' # matrix of explicit probabilities (two-dimensional test of independence)
#' gen_chisq.test(100, P0 = matrix(c(.25, .25, .25, .25), 2, 2),
#'                      P = matrix(c(.6, .2, .1, .1),2,2))
#'
gen_chisq.test <- function(n, P, ...) {
	tab <- as.vector(rmultinom(1, size = n, prob = as.vector(P)))
	if(is.matrix(P))
		tab <- matrix(tab, nrow=nrow(P), ncol=ncol(P))
	tab
}

#' Generate sample data for variance test
#'
#' Generates one or or more sets of continuous data group-level data
#' to perform a variance test.
#'
#' @param n sample size per group, assumed equal across groups
#' @param sds a vector of standard deviations to use for each group
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param ... additional arguments (not used)
#'
#' @seealso \code{\link{p_var.test}}
#' @return a data.frame with the variables 'DV' and (potentially) 'group' when
#'   constructing multi-sample data
#' @export
#' @examples
#'
#' # one sample
#' gen_var.test(100, sds=10)
#'
#' # three sample
#' gen_var.test(100, sds=c(10, 9, 11))
#'
gen_var.test <- function(n, sds, n.ratios = rep(1, length(sds)), ...){
	dat <- if(length(sds) == 1){
		data.frame(DV=rnorm(n, sd=sds))
	} else {
		dv <- sapply(1:length(sds), \(i){
			rnorm(n * n.ratios[i], sd=sds[i])
		})
		group <- rep(paste0('G', 1:length(sds)), times=n*n.ratios)
		data.frame(DV=as.vector(dv), group=group)
	}
	dat
}

#' Generate data samples for global linear regression model simulation
#'
#' Data generation for linear regression model using fixed or random
#' independent variables.
#'
#' @param n sample size
#' @param R2 R-squared effect size
#' @param k number of IVs
#' @param R2_0 null hypothesis for R-squared
#' @param k.R2_0 number of IVs associated with the null hypothesis model
#' @param R2.resid residual R-squared value, typically used when comparing
#'   nested models when fit sequentially (e.g., comparing model A vs B when
#'   model involves the structure A -> B -> C)
#' @param fixed.X logical; should the IVs be considered fixed or random?
#' @param ... additional arguments (not used)
#' @return a data.frame with the dependent variable in the left-most column
#'   with the name 'y', while all other columns reflect the independent variables
#' @seealso \code{\link{p_lm}}
#' @export
#' @examples
#'
#' # 5 fixed IVs, R^2 = .1, sample size of 95
#' gen_lm(n=95, R2=.1, k=5)
#'
gen_lm <- function(n, k, R2, R2_0 = 0, k.R2_0 = 0,
				   R2.resid=1-R2, fixed.X=TRUE, ...){
	if(!fixed.X){
		X <- matrix(rnorm(k*n), n, k)
	} else {
		lst <- vector('list', k)
		for(i in 1:k) lst[[i]] <- 0:1
		x <- expand.grid(lst)
		X <- if(nrow(x) < n)
			x[rep(1:nrow(x), each=ceiling(n/nrow(x))), ]
		else
			x[floor(seq(1, nrow(x), length.out=n)),]
		X <- X[1:n, ]
		X <- scale(X)
	}
	colnames(X) <- paste0('X', 1:k)
	R2s <- R2 - R2_0
	betas <- c(sqrt(R2s), sqrt(R2_0), rep(0, k-2))
	y <- colSums(betas * t(X)) + rnorm(n, 0, sqrt(R2.resid))
	data.frame(y, X)
}

#' Generate sample data for one-way ANOVA simulations
#'
#' Generates continuous multi-sample data to be analysed by
#' a one-way ANOVA. Samples are generated given conditional
#' observations are normally distributed and have
#' have equal variance by default, however these may be modified.
#'
#' @param n sample size per group
#' @param k number of groups
#' @param f Cohen's f effect size
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param means (optional) vector of means. When specified the input \code{f} is ignored
#' @param sds (optional) vector of SDs. When specified the input \code{f} is ignored
#' @param ... additional arguments (not used)
#' @return returns a data.frame with the variables 'DV' and 'group'
#' @seealso \code{\link{p_anova.test}}
#' @examples
#'
#' # n=50 in 3 groups, "medium" effect size
#' gen_anova.test(50, k=3, f=.25)
#'
#' # explicit means/sds
#' gen_anova.test(50, 3, means=c(0,0,1), sds=c(1,2,1))
#'
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

#' Generates three-variable mediation analysis data
#'
#' Generates continuous X and Y data and continuous or dichotomous X data
#' in a three variable mediation model, assuming Gaussian errors.
#'
#' @param n total sample size unless \code{dichotomous.X = TRUE}, in which the
#'   value represents the size per group
#' @param a regression coefficient for the path X -> M
#' @param b regression coefficient for the path M -> Y
#' @param cprime partial regression coefficient for the path X -> Y
#' @param dichotomous.X logical; should the X variable be generated as though it
#'  were dichotomous? If TRUE then \code{n} represents the sample size per group
#' @return a data.frame with the variables X, M, and Y for the independent variable,
#'  mediator, and dependent variable, respectively
#' @seealso \code{\link{p_mediate}}
#' @examples
#'
#' gen_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39)
#' gen_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, dichotomous.X=TRUE)
#'
#' @export
gen_mediation <- function(n, a, b, cprime, dichotomous.X=FALSE, ...){
	if(dichotomous.X){
		X <- rep(0:1, each=n)
		n <- n*2
	} else X <- rnorm(n)
	M <- a*X + rnorm(n)
	Y <- b*M + cprime*X + rnorm(n)
	dat <- data.frame(X, Y, M)
	dat
}

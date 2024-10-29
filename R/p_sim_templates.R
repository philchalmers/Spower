#' Independent/paired samples t-test simulation and p-value
#'
#' Generates one or two sets of continuous data group-level data
#' according to Cohen's effect size 'd', and return a p-value.
#' The data and associated t-test
#' assume that the conditional observations are normally distributed and have
#' have equal variance by default, however these may be modified.
#'
#' @param n sample size per group, assumed equal across groups
#' @param d Cohen's standardized effect size \code{d}
#' @param mu population mean to test against
#' @param type type of t-test to use; can be \code{'two.sample'},
#'   \code{'one.sample'}, or \code{'paired'}
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param var.equal logical; use the classical or Welch corrected t-test?
#' @param n2_n1 allocation ratio reflecting the same size ratio.
#'   Default of 1 sets the groups to be the same size. Only applicable
#'   when \code{type = 'two.sample'}
#' @param raw_info (optional) list of mean and SD inputs for each group,
#'   each specified as a vector. When specified the input \code{d} is ignored
#'
#'
#' @return a single p-value
#' @examples
#'
#' # sample size of 50 per group, "medium" effect size
#' p_t.test(n=50, d=0.5)
#'
#' # second group 2x as large as the first group
#' p_t.test(n=50, d=0.5, n2_n1 = 2)
#'
#' # paired and one-sample tests
#' p_t.test(n=50, d=0.5, type = 'paired')
#' p_t.test(n=50, d=0.5, type = 'one.sample')
#'
#' if(FALSE){
#'   # compare simulated results to pwr package
#'
#'   pwr::pwr.t.test(d=0.2, n=60, sig.level=0.10,
#'              type="one.sample", alternative="two.sided")
#'   Spower(n=60, d=0.2, type = 'one.sample', two.tailed=TRUE,
#'      sim=p_t.test, sig.level=.10)
#'
#'   pwr::pwr.t.test(d=0.3, power=0.75, type="two.sample",
#'                   alternative="greater")
#'   Spower(n=NA, d=0.3, type='two.sample', two.tailed=FALSE,
#'          sim=p_t.test, power=0.75, interval=c(10,200))
#'
#' }
#'
#' @export
p_t.test <- function(n, d, mu = 0,
					 type = c('two.sample', 'one.sample', 'paired'),
					 n2_n1 = 1, two.tailed = TRUE, var.equal = TRUE,
					 raw_info = list(means=NA, sds=NA)) {
	type <- match.arg(type)
	n.each <- n * n2_n1
	stopifnot(all.equal(n.each, as.integer(n.each)))
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
	p <- if(type == 'paired'){
		if(n2_n1 != 1) stop('n2_n1 must equal 1 for paired t-tests')
		t.test(group1, group2, mu=mu, paired=TRUE)$p.value
	} else if(type == 'two.sample'){
		t.test(DV ~ group, dat, var.equal=var.equal, mu=mu)$p.value
	} else if(type == 'one.sample') {
		dv <- if(!all(is.na(raw_info$means)))
			with(raw_info, rnorm(n, mean=means, sd=sds)) else rnorm(n, mean=d)
		t.test(dv, mu=mu)$p.value
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' Correlation simulation and p-value
#'
#' Generates correlated X-Y data and returns a p-value to assess the null
#' of no correlation in the population. The X-Y data are generated
#' assuming a multivariate normal distribution.
#'
#' @param n sample size
#' @param r correlation
#' @param rho population coefficient to test against using
#'   \code{\link[car]{linearHypothesis}}. Defaults to 0
#' @param method method to use to compute the correlation
#'   (see \code{\link{cor.test}}). Only used when \code{rho = 0}
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @importFrom car linearHypothesis
#' @return a single p-value
#' @examples
#'
#' # 50 observations, .5 correlation
#' p_r(50, r=.5)
#' p_r(50, r=.5, method = 'spearman')
#'
#' # test against constant other than rho = .6
#' p_r(50, .5, rho=.60)
#'
#'
#' if(FALSE){
#'     # compare simulated results to pwr package
#'
#'     pwr::pwr.r.test(r=0.3, n=50)
#'     Spower(n=50, r=0.3, sim=p_r)
#'
#'     pwr::pwr.r.test(r=0.3, power=0.80)
#'     Spower(n=NA, r=0.3, sim=p_r, power=.80, interval=c(10, 200))
#'
#'     pwr::pwr.r.test(r=0.1, power=0.80)
#'     Spower(n=NA, r=0.1, sim=p_r, power=.80, interval=c(200, 1000))
#'
#' }
#'
#' @export
p_r <- function(n, r, rho = 0, method = 'pearson', two.tailed = TRUE) {
	dat <- SimDesign::rmvnorm(n, sigma = matrix(c(1,r,r,1), 2, 2))
	colnames(dat) <- c('x', 'y')
	if(rho != 0){
		mod <- lm(y ~ x, as.data.frame(dat))
		out <- car::linearHypothesis(mod, "x", rho)
		p <- out$`Pr(`[2]
	} else {
		p <- cor.test(~ x + y, dat, method=method)$p.value
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' Proportion test simulation and p-value
#'
#' Generates single and multi-sample data
#' for proportion tests and return a p-value. Uses \code{\link{binom.test}}
#' for one-sample applications and \code{\link{prop.test}} otherwise.
#'
#' @param n sample size per group
#' @param prop sample probability/proportions of success.
#'   If a vector with two-values or more elements are supplied then
#'   a multi-samples test will be used
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param pi probability of success to test against (default is .5). Ignored
#'   for two-sample tests
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @return a single p-value
#' @examples
#'
#' # 50 observations, test against pi = .5
#' p_prop.test(50, prop=.65)
#'
#' # two sample test
#' p_prop.test(50, prop=c(.5, .65))
#'
#' # two sample test, unequal ns
#' p_prop.test(50, prop=c(.5, .65))
#'
#' # three sample test, group2 twice as large as others
#' p_prop.test(50, prop=c(.5, .65, .7), n.ratios=c(1,2,1))
#'
#' if(FALSE){
#'     # compare simulated results to pwr package
#'
#'     h <- pwr::ES.h(0.5, 0.4)
#'     pwr::pwr.p.test(h=h, n=60)
#'
#'     Spower(n=60, prop=c(.5, .4), sim=p_prop.test)
#'
#' }
#'
#' @export
p_prop.test <- function(n, prop, pi = .5,
						n.ratios = rep(1, length(prop)),
						two.tailed = TRUE) {
	stopifnot(length(n) == 1)
	n.each <- n * n.ratios
	stopifnot(all.equal(n.each, as.integer(n.each)))
	p <- if(length(prop) > 1){
		draws <- sapply(1:length(prop), \(i){
			vals <- rbinom(n * n.ratios[i], 1, prob=prop[i])
			c(sum(vals), length(vals))
		})
		A <- draws[1,]
		B <- draws[2,]
		prop.test(A, B)$p.value
	} else {
		dat <- rbinom(n, 1, prob = prop)
		binom.test(sum(dat), n=n, p=pi)$p.value
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' One-way ANOVA simulation and p-value
#'
#' Generates continuous multi-sample data to be analysed by
#' a one-way ANOVA, and return a p-value.
#' Uses the function \code{\link{one.way}} to perform the analyses.
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
#' @param raw_info (optional) list of mean and SD inputs for each group,
#'   each specified as a vector. When specified the input \code{f} is ignored
#' @return a single p-value
#' @examples
#'
#' # n=50 in 3 groups, "medium" effect size
#' p_anova.test(50, k=3, f=.25)
#'
#' # explicit means/sds
#' p_anova.test(50, 3,
#'             raw_info=list(means=c(0,0,1), sds=c(1,2,1)))
#'
#' if(FALSE){
#'   # compare simulated results to pwr package
#'   pwr::pwr.anova.test(f=0.28, k=4, n=20)
#'   Spower(n=20, k=4, f=.28, sim=p_anova.test)
#' }
#'
#' @export
p_anova.test <- function(n, k, f,
						 n.ratios = rep(1, k),
						 two.tailed = TRUE, var.equal = TRUE,
						 raw_info = list(means=NA, sds=NA)) {
	stopifnot(length(n) == 1)
	stopifnot(length(n.ratios) == k)
	group <- rep(factor(1:k), times = n*n.ratios)
	n.each <- n*n.ratios
	stopifnot(all.equal(n.each, as.integer(n.each)))
	if(!all(is.na(raw_info$means))){
		dv <- sapply(1:k, \(i)
					  with(raw_info, rnorm(n*n.ratios[i],
					  					 mean=means[i], sd=sds[i])))
		df <- data.frame(group=group, dv = as.numeric(dv))
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
		df <- data.frame(dv=dv, group=group)
	}
	p <- oneway.test(dv ~ group, data=df, var.equal=var.equal)$p.value
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' Proportion test simulation and p-value
#'
#' Generates single and multi-sample data
#' for proportion tests and return a p-value. Uses \code{\link{binom.test}}
#' for one-sample applications and \code{\link{prop.test}} otherwise.
#'
#' @param n sample size per group
#' @param prop sample probability/proportions of success.
#'   If a vector with two-values or more elements are supplied then
#'   a multi-samples test will be used
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param pi probability of success to test against (default is .5). Ignored
#'   for two-sample tests
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @return a single p-value
#' @examples
#'
#' p_chisq.test(100, w=.2, df=3)
#'
#' p_chisq.test(100, raw_info = list(P0 = c(.25, .25, .25, .25),
#'                                   P = c(.6, .2, .1, .1)))
#'
#' # works, but logic seems odd
#' p_chisq.test(100, raw_info = list(P0 = matrix(c(.25, .25, .25, .25), 2, 2),
#'                                   P = matrix(c(.6, .2, .1, .1),2,2)))
#'
#' if(FALSE){
#'     # compare simulated results to pwr package
#'
#'     P0 <- c(1/3, 1/3, 1/3)
#'     P <- c(.5, .25, .25)
#'     w <- pwr::ES.w1(P0, P)
#'     df <- 3-1
#'     pwr::pwr.chisq.test(w=w, df=df, N=100, sig.level=0.05)
#'
#'     Spower(n=100, w=w, df=df, sim=p_chisq.test)
#'     Spower(n=100, sim=p_chisq.test, raw_info=list(P0=P0, P=P))
#'
#'
#'
#' }
#'
#' @export
p_chisq.test <- function(n, w, df,
						 raw_info = list(P0 = NA, P = NA)) {
	stopifnot(length(n) == 1)
	p <- if(!missing(w)){
		stopifnot(length(w) == 1)
		stopifnot(length(df) == 1)
		w2 <- w^2
		X2 <- n*w2
		X2_n <- X2/n
		p0 <- 1/(df+1)
		P <- abs(p0 - sqrt(X2_n * p0))
		#n * (P - p0)^2 / p0
		tab <- as.vector(rmultinom(1, size = n, prob = rep(P, df+1)))
		p <- chisq.test(tab, p=rep(p0, df+1))$p.value
	} else {
		tab <- as.vector(with(raw_info, rmultinom(1, size = n, prob = P)))
		if(is.matrix(raw_info$P))
			tab <- with(raw_info, matrix(tab, nrow=nrow(P), ncol=ncol(P)))
		p <- chisq.test(tab, p=raw_info$P0)$p.value
	}
	p
}

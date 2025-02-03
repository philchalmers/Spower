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
#' @param r (optional) instead of specifying \code{d} specify
#'   a point-biserial correlation. Internally this is transformed
#'   into a suitable \code{d} value for the power computations
#' @param type type of t-test to use; can be \code{'two.sample'},
#'   \code{'one.sample'}, or \code{'paired'}
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param var.equal logical; use the classical or Welch corrected t-test?
#' @param n2_n1 allocation ratio reflecting the same size ratio.
#'   Default of 1 sets the groups to be the same size. Only applicable
#'   when \code{type = 'two.sample'}
#' @param raw_info (optional) list of mean and SD inputs for each group,
#'   each specified as a vector. When specified the input \code{d} is ignored
#' @param gen_fun function used to generate the required two-sample data.
#'   Object returned must be a \code{data.frame} with the columns
#'   \code{"DV"} and \code{"group"}. Default uses \code{\link{gen_t.test}}
#'   to generate conditionally Gaussian distributed samples
#'
#'   User defined version of this function must, at minimum, accept the arguments
#'   \code{gen_fun(n, n2_n1, d, r, type, raw_info, ...)}, even if they are not used
#'   directly
#' @param ... additional arguments to be passed to \code{gen_fun}. Not necessary
#'   unless a customized \code{gen_fun} was defined
#'
#'
#' @return a single p-value
#' @examples
#'
#' # sample size of 50 per group, "medium" effect size
#' p_t.test(n=50, d=0.5)
#'
#' # point-biserial correlation effect size
#' p_t.test(n=50, r=.3)
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
#'   Spower(p_t.test, n=60, d=0.2, type = 'one.sample', two.tailed=TRUE,
#'          sig.level=.10)
#'
#'   pwr::pwr.t.test(d=0.3, power=0.80, type="two.sample",
#'                   alternative="greater")
#'   Spower(p_t.test, n=NA, d=0.3, type='two.sample', two.tailed=FALSE,
#'          power=0.80, interval=c(10,200))
#'
#' }
#'
#' ###### Custom data generation function
#'
#' # Generate data such that:
#' #   - group 1 is from a negatively distribution (reversed X2(10)),
#' #   - group 2 is from a positively skewed distribution (X2(5))
#' #   - groups have equal variance, but differ by d = 0.5
#'
#' args(gen_t.test)   ## use these arguments as a basis
#'
#' my.gen_fun <- function(n, d, n2_n1, r, type, df1, df2, ...){
#'  	 group1 <- -1 * rchisq(n, df=df1)
#' 	     group2 <- rchisq(n, df=df2)
#' 	     # scale groups first given moments of the chi-square distribution,
#' 	     #   then add std mean difference
#' 	     group1 <- ((group1 + df1) / sqrt(2*df1))
#' 	     group2 <- ((group2 - df2) / sqrt(2*df2)) + d
#' 	     dat <- data.frame(DV=c(group1, group2),
#' 	        			   group=gl(2, n, labels=c('G1', 'G2')))
#' 	     dat
#' }
#' df <- my.gen_fun(n=10000, d=.5, df1=10, df2=5)
#' with(df, tapply(DV, group, mean))
#' with(df, tapply(DV, group, sd))
#'
#' library(ggplot2)
#' ggplot(df, aes(group, DV, fill=group)) + geom_violin()
#'
#' p_t.test(n=100, d=0.5, gen_fun=my.gen_fun, df1=10, df2=5)
#'
#' if(FALSE){
#'
#'   # power given Gaussian distributions
#'   Spower(p_t.test, n=100, d=0.5)
#'
#'   # estimate power given the customized data generating function
#'   Spower(p_t.test, n=100, d=0.5, gen_fun=my.gen_fun, df1=10, df2=5)
#'
#'   # evaluate Type I error rate to see if liberal/conservative given
#'   # assumption violations (should be close to alpha/sig.level)
#'   Spower(p_t.test, n=100, d=0, gen_fun=my.gen_fun, df1=10, df2=5)
#'
#' }
#'
#' @export
p_t.test <- function(n, d, mu = 0, r = NULL,
					 type = c('two.sample', 'one.sample', 'paired'),
					 n2_n1 = 1, two.tailed = TRUE, var.equal = TRUE,
					 raw_info = list(means=NA, sds=NA),
					 gen_fun=gen_t.test, ...) {
	type <- match.arg(type)
	if(!missing(d) && !is.null(r))
		stop('Please use either d or r')
	if(!is.null(r)){
		type <- 'two.sample'
		stopifnot(var.equal)
	}
	dat <- gen_fun(n=n, n2_n1=n2_n1, d=d, r=r, type=type,
				   raw_info=raw_info, ...)
	p <- if(type == 'paired'){
		if(n2_n1 != 1) stop('n2_n1 must equal 1 for paired t-tests')
		lvls <- levels(dat$group)
		group1 <- with(dat, DV[group == lvls[1]])
		group2 <- with(dat, DV[group == lvls[2]])
		t.test(group1, group2, mu=mu, paired=TRUE)$p.value
	} else if(type == 'two.sample'){
		t.test(DV ~ group, dat, var.equal=var.equal, mu=mu)$p.value
	} else if(type == 'one.sample') {
		t.test(dat$DV, mu=mu)$p.value
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' Generate sample data for t.test procedures
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

#' Correlation simulation and p-value
#'
#' Generates correlated X-Y data and returns a p-value to assess the null
#' of no correlation in the population. The X-Y data are generated
#' assuming a multivariate normal distribution.
#'
#' @param n sample size
#' @param r correlation
#' @param rho population coefficient to test against. Uses the
#'   Fisher's z-transformation approximation when non-zero
#' @param method method to use to compute the correlation
#'   (see \code{\link{cor.test}}). Only used when \code{rho = 0}
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
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
#'     Spower(p_r, n=50, r=0.3)
#'
#'     pwr::pwr.r.test(r=0.3, power=0.80)
#'     Spower(p_r, n=NA, r=0.3, power=.80, interval=c(10, 200))
#'
#'     pwr::pwr.r.test(r=0.1, power=0.80)
#'     Spower(p_r, n=NA, r=0.1, power=.80, interval=c(200, 1000))
#'
#' }
#'
#' @export
p_r <- function(n, r, rho = 0, method = 'pearson', two.tailed = TRUE) {
	dat <- SimDesign::rmvnorm(n, sigma = matrix(c(1,r,r,1), 2, 2))
	colnames(dat) <- c('x', 'y')
	if(rho != 0){
		out <- cor.test(~ x + y, dat, method=method)
		z <- with(out, 1/2 * log((1+estimate)/(1-estimate)))
		se <- 1 / sqrt(n-3)
		z0 <- 1/2 * log((1+rho)/(1-rho))
		t <- (z - z0) / se
		p <- pnorm(abs(t), lower.tail=FALSE)*2
	} else {
		p <- cor.test(~ x + y, dat, method=method)$p.value
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' Polychoric and polyserial simulation and p-value
#'
#' Generates correlated X-Y data and returns a p-value to assess the null
#' of no correlation in the population. The X-Y data are generated
#' assuming a multivariate normal distribution and subsequently discretized
#' for one or both of the variables.
#'
#' @param n sample size
#' @param r correlation prior to the discretization (recovered via the
#'   polyserial/polychoric estimates)
#' @param rho population coefficient to test against
#' @param tauX intercept parameters used for discretizing the X variable
#' @param tauY intercept parameters used for discretizing the Y variable. If
#'   missing a polyserial correlation will be estimated, otherwise a
#'   tetrachoric/polychoric correlation will be estimated
#' @param ML logical; use maximum-likelihood estimation?
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param score logical; should the SE be based at the null hypothesis (score test)
#'   or the ML estimate (Wald test)? The former is the canonical form for
#'   a priori power analyses though requires twice as many computations as the
#'   Wald test approach
#'
#' @return a single p-value
#' @export
#' @examples
#'
#' # 100 observations, .5 correlation, tetrachoric estimate
#' p_r.cat(100, r=.5, tauX=0, tauY=1)
#'
#' # Wald test
#' p_r.cat(100, r=.5, tauX=0, tauY=1, score=FALSE)
#'
#' # polyserial estimate (Y continuous)
#' p_r.cat(50, r=.5, tauX=0)
#'
p_r.cat <- function(n, r, tauX, rho=0, tauY = NULL,
					ML=TRUE, two.tailed=TRUE, score=FALSE){
	continuous.Y <- is.null(tauY)
	dat <- SimDesign::rmvnorm(n,
							  sigma = matrix(c(1,r,r,1), 2, 2))
	datcut <- matrix(0, n, 2)
	for(i in length(tauX):1)
		datcut[dat[,1] > tauX[i], 1] <- i
	if(!continuous.Y){
		for(i in length(tauY):1)
			datcut[dat[,2] > tauY[i], 2] <- i
	} else datcut[,2] <- dat[,2]
	datcut <- as.data.frame(datcut)
	colnames(datcut) <- c('x', 'y')
	# could generate r=0 data to get SE_0 instead for proper score test
	out <- if(continuous.Y){
		with(datcut, polycor::polyserial(y, x, ML=ML, std.err=TRUE))
	} else {
		with(datcut, polycor::polychor(y, x, ML=ML, std.err=TRUE))
	}
	est <- out$rho
	vcov <- out$var
	if(score > 1) return(sqrt(vcov[1,1]))
	SE <- if(score == 1){
		p_r.cat(n, r=rho, tauX=tauX, rho=rho, tauY=tauY,
				ML=ML, two.tailed=two.tailed, score=2)
	} else sqrt(vcov[1,1])
	z <- (est - rho) / SE
	p <- pnorm(abs(z), lower.tail=FALSE)
	p <- ifelse(two.tailed, p*2, p)
	p
}

if(FALSE){
	p_r.cat(100, .3, tauX=0, tauY=1)

	p_r.cat(100, .3, tauX=0)

	Spower(p_r.cat, n=100, r=.3, tauX=0, tauY=1, parallel=TRUE)



}



#' Proportion test simulation and p-value
#'
#' Generates single and multi-sample data
#' for proportion tests and return a p-value. Uses \code{\link{binom.test}}
#' for one-sample applications and \code{\link{prop.test}} otherwise.
#'
#' @param n sample size per group
#' @param h Cohen's h effect size; only supported for one-sample analysis.
#'
#'   Note that it's important to specify the null
#'   value \code{pi} when supplying this effect size as the power
#'   changes depending on these specific values (see example below).
#' @param prop sample probability/proportions of success.
#'   If a vector with two-values or more elements are supplied then
#'   a multi-samples test will be used. Matrices are also supported
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param pi probability of success to test against (default is .5). Ignored
#'   for two-sample tests
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param exact logical; use fisher's exact test via \code{\link{fisher.test}}?
#'   Use of this flag requires that \code{prop} was specified as a matrix
#' @param correct logical; use Yates' continuity correction?
#' @return a single p-value
#' @examples
#'
#' # one sample, 50 observations, tested against pi = .5 by default
#' p_prop.test(50, prop=.65)
#'
#' # specified using h and pi
#' h <- pwr::ES.h(.65, .4)
#' p_prop.test(50, h=h, pi=.4)
#' p_prop.test(50, h=-h, pi=.65)
#'
#' # two-sample test
#' p_prop.test(50, prop=c(.5, .65))
#'
#' # two-sample test, unequal ns
#' p_prop.test(50, prop=c(.5, .65), n.ratios = c(1,2))
#'
#' # three-sample test, group2 twice as large as others
#' p_prop.test(50, prop=c(.5, .65, .7), n.ratios=c(1,2,1))
#'
#' # Fisher exact test
#' p_prop.test(50, prop=matrix(c(.5, .65, .7, .5), 2, 2))
#'
#' if(FALSE){
#'     # compare simulated results to pwr package
#'
#'     # one-sample tests
#'     (h <- pwr::ES.h(0.5, 0.4))
#'     pwr::pwr.p.test(h=h, n=60)
#'
#'     # uses binom.test (need to specify null location as this matters!)
#'     Spower(p_prop.test, n=60, h=h, pi=.4)
#'     Spower(p_prop.test, n=60, prop=.5, pi=.4)
#'
#'     # compare with switched null
#'     Spower(p_prop.test, n=60, h=h, pi=.5)
#'     Spower(p_prop.test, n=60, prop=.4, pi=.5)
#'
#'     # two-sample test, one-tailed
#'     (h <- pwr::ES.h(0.67, 0.5))
#'     pwr::pwr.2p.test(h=h, n=80, alternative="greater")
#'     Spower(p_prop.test, n=80, prop=c(.67, .5), two.tailed=FALSE, correct=FALSE)
#'
#'     # same as above, but with continuity correction (default)
#'     Spower(p_prop.test, n=80, prop=c(.67, .5), two.tailed=FALSE)
#'
#'     # three-sample joint test, equal n's
#'     Spower(p_prop.test, n=50, prop=c(.6,.4,.7))
#'
#' }
#'
#' @export
p_prop.test <- function(n, h, prop, pi = .5,
						n.ratios = rep(1, length(prop)),
						two.tailed = TRUE, correct=TRUE, exact=FALSE) {
	root.h <- function(p1, p2, h)
		(2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))) - h
	stopifnot(length(n) == 1)
	if(!missing(h)){
		n.ratios <- 1
		int <- if(h > 0) c(pi, 1) else c(0, pi)
		root <- uniroot(root.h, interval=int,
						p2=ifelse(length(n.ratios) == 2, prop, pi), h=h)$root
		prop <- root
		# pwr::ES.h(prop, pi) == h
	}
	n.each <- n * n.ratios
	stopifnot(all.equal(n.each, as.integer(n.each)))
	p <- if(length(prop) > 1){
		draws <- sapply(1:length(prop), \(i){
			vals <- rbinom(n * n.ratios[i], 1, prob=prop[i])
			c(sum(vals), length(vals))
		})
		A <- draws[1,]
		B <- draws[2,]
		if(exact){
			stopifnot(is.matrix(prop))
			A <- matrix(A, nrow(prop), ncol(prop))
			fisher.test(A)$p.value
		} else prop.test(A, B, correct=correct)$p.value
	} else {
		dat <- rbinom(n, 1, prob = prop)
		binom.test(sum(dat), n=n, p=pi)$p.value
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' McNemar test simulation and p-value
#'
#' Generates two-dimensional sample data for McNemar test and
#' return a p-value. Uses \code{\link{mcnemar.test}}.
#'
#' @param n total sample size
#' @param prop two-dimensional matrix of proportions/probabilities
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
# @param exact logical; use fisher's exact test via \code{\link{fisher.test}}?
#' @param correct logical; use continuity correction? Only applicable for
#'   2x2 tables
#' @return a single p-value
#' @examples
#'
#' # from ?mcnemar.test
#' Performance <- matrix(c(794, 86, 150, 570),
#' 		   nrow = 2,
#' 		   dimnames = list("1st Survey" = c("Approve", "Disapprove"),
#' 		               "2nd Survey" = c("Approve", "Disapprove")))
#' (prop <- prop.table(Performance))
#'
#' # one sample + test and resulting p-value
#' p_mcnemar.test(n=sum(Performance), prop=prop)
#'
#' \dontrun{
#'
#' # post-hoc power (not recommended)
#' Spower(p_mcnemar.test, n=sum(Performance), prop=prop)
#'
#' }
#'
#' @export
p_mcnemar.test <- function(n, prop,
						   two.tailed = TRUE, correct=TRUE) {
	draws <- rmultinom(1, n, prob = as.numeric(prop))
	p <- mcnemar.test(matrix(draws, nrow(prop), ncol(prop)),
					  correct=correct)$p.value
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' One-way ANOVA simulation and p-value
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
#'   Spower(p_anova.test, n=20, k=4, f=.28)
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

#' Chi-squared test simulation and p-value
#'
#' Generates multinomial data suitable for analysis with
#' \code{\link{chisq.test}}.
#'
#' @param n sample size per group
#' @param w Cohen's w effect size
#' @param df degrees of freedom
#' @param correct logical; apply continuity correction?
#' @param raw_info list of raw information to generate multinomial data
#'   under specific nulls and power configurations
#'
#' @return a single p-value
#' @examples
#'
#' # effect size w + df
#' p_chisq.test(100, w=.2, df=3)
#'
#' # vector of explicit probabilities
#' p_chisq.test(100, raw_info = list(P0 = c(.25, .25, .25, .25),
#'                                   P = c(.6, .2, .1, .1)))
#'
#' # matrix of explicit probabilities
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
#'     Spower(p_chisq.test, n=100, w=w, df=df)
#'     Spower(p_chisq.test, n=100, raw_info=list(P0=P0, P=P))
#'
#'     # slightly differ (latter more conservative)
#'     pwr::pwr.chisq.test(w=w, df=df, power=.8, sig.level=0.05)
#'     Spower(p_chisq.test, n=NA, w=w, df=df,
#'            power=.80, interval=c(50, 200))
#'
#'     # Spower more conservative even with large N
#'     pwr::pwr.chisq.test(w=.1, df=df, power=.95, sig.level=0.05)
#'     out <- Spower(p_chisq.test, n=NA, w=.1, df=df,
#'                  power=.95, interval=c(1000, 2000))
#'     summary(out)
#'
#' }
#'
#' @export
p_chisq.test <- function(n, w, df,
						 correct = TRUE,
						 raw_info = list(P0 = NA, P = NA)) {
	stopifnot(length(n) == 1)
	p <- if(!missing(w)){
		stopifnot(length(w) == 1)
		stopifnot(length(df) == 1)
		w2 <- w^2
		p0 <- rep(1/(df+1), df+1)
		fn <- function(p1, p0, df, w2){
			ps <- c(p1, rep((1 - p1)/ df, df))
			(sum((ps - p0)^2 / p0) - w2)^2
		}
		# strange that optimize(fn, c(0,1)) gives right w2 but wrong p?
		opt <- optimize(fn, c(p0[1],1), p0=p0, df=df, w2=w2)
		P <- with(opt, c(minimum, rep((1 - minimum)/df, df)))
		# sum((P - p0)^2 / p0) # == w2
		tab <- as.vector(rmultinom(1, size = n, prob = P))
		p <- chisq.test(tab, correct=correct, p=p0)$p.value
	} else {
		tab <- as.vector(with(raw_info, rmultinom(1, size = n, prob = P)))
		if(is.matrix(raw_info$P))
			tab <- with(raw_info, matrix(tab, nrow=nrow(P), ncol=ncol(P)))
		p <- chisq.test(tab, correct=correct, p=raw_info$P0)$p.value
	}
	p
}

#' Variance test simulation and p-value
#'
#' Generates one or or more sets of continuous data group-level data
#' to perform a variance test, and return a p-value. When two-samples
#' are investigated the \code{\link{var.test}} function will be used,
#' otherwise functions from the \code{EnvStats} package will be used.
#'
#'
#' @param n sample size per group, assumed equal across groups
#' @param sds a vector of standard deviations to use for each group
#' @param ratio hypothesized ratio of the population variance (only
#'   used in two-sample case)
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param test type of test to use in multi-sample applications.
#'   Can be either \code{'Levene'} (default) or \code{'Bartlett'}
#' @param sigma standard deviation value to test against in one-sample test
#' @param correct logical; use correction when \code{test = 'Bartlett'}?
#'
#' @return a single p-value
#' @export
#' @importFrom EnvStats varTest varGroupTest
#' @examples
#'
#' # one sample
#' p_var.test(100, sds=10, sigma=9)
#'
#' # three sample
#' p_var.test(100, sds=c(10, 9, 11))
#'
#' if(FALSE){
#'   # power to detect three-group variance differences
#'   Spower(p_var.test, n=100, sds=c(10,9,11))
#'
#'   # sample size per group to achieve 80% power
#'   Spower(p_var.test, n=NA, sds=c(10,9,11),
#'          power=.80, interval=c(100, 1000))
#' }
#'
p_var.test <- function(n, sds, n.ratios = rep(1, length(sds)),
					   ratio = 1, sigma = 1, two.tailed = TRUE,
					   test = 'Levene', correct = TRUE){
	p <- if(length(sds) == 1){
		dv <- rnorm(n, sd=sds)
		out <- EnvStats::varTest(dv, sigma.squared = sigma^2)
		out$p.value
	} else {
		dv <- sapply(1:length(sds), \(i){
			rnorm(n * n.ratios[i], sd=sds[i])
		}) |> as.vector()
		group <- rep(paste0('G', 1:length(sds)), times=n*n.ratios)
		out <- if(length(sds) == 2) var.test(dv ~ group)
			else EnvStats::varGroupTest(dv ~ group, test=test,
									  correct=correct)
		out$p.value
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' Wilcox test simulation and p-value
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

#' Linear regression model and p-value
#'
#' p-values associated with linear regression model using fixed or random
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
#' @export
#' @examples
#'
#' # 5 fixed IVs, R^2 = .1, sample size of 95
#' p_lm(n=95, R2=.1, k=5)
#'
p_lm <- function(n, R2, k, R2_0 = 0, k.R2_0 = 0, R2.resid=1-R2, fixed.X=TRUE){
	stopifnot(R2 > R2_0)
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
	if(k.R2_0 == 0)
		stopifnot(k >= 2)
	R2s <- R2 - R2_0
	betas <- c(sqrt(R2s), sqrt(R2_0), rep(0, k-2))
	y <- colSums(betas * t(X)) + rnorm(n, 0, sqrt(R2.resid))
	df <- data.frame(y, X)
	if(k.R2_0 > 0)
		df2 <- df[,c(1,3:(k.R2_0 +2))]
	mod1 <- lm(y ~ ., df)
	if(!fixed.X && k.R2_0 == 0 && R2_0 != 0){
		stop('Random X with non-zero R2_0 not currently supported', call.=FALSE)
	} else {
		mod0 <- if(R2_0 == 0) lm(y ~ 1, df)
		  else lm(y ~ ., df2)
	}
	p <- anova(mod0, mod1)[2, "Pr(>F)"]
	p
}

if(FALSE){

	# Example 7.3b
	# G*power gives 0.3464 (broken)
	Spower(p_lm, n=100, R2=.4, R2_0 = .3, k=5, fixed.X=FALSE)

	# Example 7.3c
	# G*power gives N=153 (broken)
	Spower(p_lm, n=NA, R2=.05, R2_0 = .2, k=5, fixed.X=FALSE, power=.9,
		   interval=c(50,300))


}

# example should include Poisson distribution to match
p_glm <- function(n){

}

#' Compare two correlations and p-value
#'
#' Function utilizes \code{\link[cocor]{cocor}} to perform correlation
#' comparison for independent, overlapping, and non-overlapping designs.
#'
#' @param n sample size
#' @param r.ab1 correlation between variable A and B in sample 1
#' @param r.ab2 correlation between variable A and B in sample 2
#' @param r.ac1 same pattern as \code{r.ab1}
#' @param r.ac2 same pattern as \code{r.ab2}
#' @param r.bc1 ...
#' @param r.bc2 ...
#' @param r.ad1 ...
#' @param r.ad2 ...
#' @param r.bd1 ...
#' @param r.bd2 ...
#' @param r.cd1 ...
#' @param r.cd2 ...
#' @param n2_n1 sample size ratio
#' @param two.tailed logical; use two-tailed test?
#' @param type type of correlation design
#' @param test hypothesis method to use. Defaults to 'fisher1925'
#' @importFrom cocor cocor
#' @importFrom methods slot
#' @export
#' @examples
#'
#' # independent (same x-y pairing across groups)
#' p_2r(100, r.ab1=.5, r.ab2=.6)
#'
#' # overlap (same y, different xs)
#' p_2r(100, r.ab1=.5, r.ab2=.7,
#'           r.ac1=.3, r.ac2=.3,
#'           r.bc1=.2, r.bc2=.2, type = 'overlap')
#'
#' # nonoverlap (different ys, different xs)
#' p_2r(100, r.ab1=.5, r.ab2=.6,
#'           r.ac1=.3, r.ac2=.3,
#'           r.bc1=.2, r.bc2=.2,
#'           r.ad1=.2, r.ad2=.2,
#'           r.bd1=.4, r.bd2=.4,
#'           r.cd1=.2, r.cd2=.2,
#'           type = 'nonoverlap')
#'
#'
p_2r <- function(n, r.ab1, r.ab2, r.ac1, r.ac2, r.bc1, r.bc2,
				 r.ad1, r.ad2, r.bd1, r.bd2, r.cd1, r.cd2,
				 n2_n1 = 1, two.tailed=TRUE,
				 type = c('independent', 'overlap', 'nonoverlap'),
				 test = 'fisher1925'){
	type <- match.arg(type)
	if(type == 'independent'){
		R1 <- matrix(c(1,r.ab1, r.ab1, 1), 2, 2)
		R2 <- matrix(c(1,r.ab2, r.ab2, 1), 2, 2)
		cnms <- c('y', 'x')
	} else if(type == 'overlap'){
		R1 <- matrix(c(1,r.ab1, r.ac1,
					   r.ab1, 1, r.bc1,
					   r.ac1,r.bc1, 1), 3, 3)
		R2 <- matrix(c(1,r.ab2, r.ac2,
					   r.ab2, 1, r.bc2,
					   r.ac2, r.bc2, 1), 3, 3)
		cnms <- c('y', 'x1', 'x2')
	} else {
		R1 <- matrix(c(1,r.ab1, r.ac1, r.ad1,
					   0, 1, r.bc1, r.bd1,
					   0,0, 1, r.cd1,
					   0, 0, 0, 1), 4, 4)
		R1 <- R1 + t(R1) - diag(4)
		R2 <- matrix(c(1,r.ab2, r.ac2, r.ad2,
					   0, 1, r.bc2, r.bd2,
					   0,0, 1, r.cd2,
					   0, 0, 0, 1), 4, 4)
		R2 <- R2 + t(R2) - diag(4)
		cnms <- c('y1', 'x1', 'y2', 'x2')
	}
	df1 <- data.frame(rmvnorm(n, sigma=R1))
	df2 <- data.frame(rmvnorm(n * n2_n1, sigma=R2))
	colnames(df1) <- colnames(df2) <- cnms
	dat <- list(sample1=df1, sample2=df2)
	res <- if(type == 'independent'){
		cocor::cocor(~ y + x | y + x, dat, test=test)
	} else if(type == 'overlap'){
		cocor::cocor(~ y + x1 | y + x2, dat, test=test)
	} else {
		cocor::cocor(~ y1 + x1 | y2 + x2, dat, test=test)
	}
	pick <- methods::slot(res, test)
	p <- pick$p.value
	p <- ifelse(two.tailed, p, p*2)
	p
}


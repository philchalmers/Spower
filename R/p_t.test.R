#' p-value from independent/paired samples t-test simulation
#'
#' Generates one or two sets of continuous data group-level data
#' according to Cohen's effect size 'd', and returns a p-value.
#' The data and associated t-test
#' assume that the conditional observations are normally distributed and have
#' have equal variance by default, however these may be modified.
#'
#' @param n sample size per group, assumed equal across groups. For paired samples
#'   this corresponds to the number of pairs (hence, half the number of data points
#'   observed)
#' @param d Cohen's standardized effect size \code{d}. For the generated data this standardized
#'   mean appears in the first group (two-sample)/first time point (paired samples)
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
#' @param means (optional) vector of means for each group.
#'   When specified the input \code{d} is ignored
#' @param sds (optional) vector of SDs for each group. If not specified and \code{d}
#'  is used then these are set to a vector of 1's
#' @param conf.level confidence interval level passed
#'   to \code{\link[stats]{t.test}}
#' @param gen_fun function used to generate the required two-sample data.
#'   Object returned must be a \code{list} containing one (one-sample) or
#'   two (independent samples/paired samples) elements,
#'   both of which are \code{numeric} vectors. Default uses \code{\link{gen_t.test}}
#'   to generate conditionally Gaussian distributed samples.
#'   User defined version of this function must include the argument \code{...}
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#'
#' @seealso \code{\link{gen_t.test}}
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
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
#' # specify mean/SDs explicitly
#' p_t.test(n=50, means = c(0,1), sds = c(2,2))
#'
#' # paired and one-sample tests
#' p_t.test(n=50, d=0.5, type = 'paired') # n = number of pairs
#' p_t.test(n=50, d=0.5, type = 'one.sample')
#'
#' # return analysis object
#' p_t.test(n=50, d=0.5, return_analysis=TRUE)
#'
#' \donttest{
#'   # compare simulated results to pwr package
#'
#'   pwr::pwr.t.test(d=0.2, n=60, sig.level=0.10,
#'              type="one.sample", alternative="two.sided")
#'   p_t.test(n=60, d=0.2, type = 'one.sample', two.tailed=TRUE) |>
#'          Spower(sig.level=.10)
#'
#'   pwr::pwr.t.test(d=0.3, power=0.80, type="two.sample",
#'                   alternative="greater")
#'   p_t.test(n=NA, d=0.3, type='two.sample', two.tailed=FALSE) |>
#'          Spower(power=0.80, interval=c(10,200))
#'
#' }
#'
#'
#' ###### Custom data generation function
#'
#' # Generate data such that:
#' #   - group 1 is from a negatively distribution (reversed X2(10)),
#' #   - group 2 is from a positively skewed distribution (X2(5))
#' #   - groups have equal variance, but differ by d = 0.5
#'
#' args(gen_t.test)   ## can use these arguments as a basis, though must include ...
#'
#' # arguments df1 and df2 added; unused arguments caught within ...
#' my.gen_fun <- function(n, d, df1, df2, ...){
#'  	 group1 <- -1 * rchisq(n, df=df1)
#' 	     group2 <- rchisq(n, df=df2)
#' 	     # scale groups first given moments of the chi-square distribution,
#' 	     #   then add std mean difference
#' 	     group1 <- ((group1 + df1) / sqrt(2*df1))
#' 	     group2 <- ((group2 - df2) / sqrt(2*df2)) + d
#' 	     dat <- list(group1, group2)
#' 	     dat
#' }
#'
#' # check the sample data properties
#' dat <- my.gen_fun(n=10000, d=.5, df1=10, df2=5)
#' sapply(dat, mean)
#' sapply(dat, sd)
#'
#' p_t.test(n=100, d=0.5, gen_fun=my.gen_fun, df1=10, df2=5)
#'
#' \donttest{
#'
#'   # power given Gaussian distributions
#'   p_t.test(n=100, d=0.5) |> Spower(replications=30000)
#'
#'   # estimate power given the customized data generating function
#'   p_t.test(n=100, d=0.5, gen_fun=my.gen_fun, df1=10, df2=5) |>
#'     Spower(replications=30000)
#'
#'   # evaluate Type I error rate to see if liberal/conservative given
#'   # assumption violations (should be close to alpha/sig.level)
#'   p_t.test(n=100, d=0, gen_fun=my.gen_fun, df1=10, df2=5) |>
#'     Spower(replications=30000)
#'
#' }
#'
#' @export
p_t.test <- function(n, d, mu = 0, r = NULL, type = 'two.sample',
					 n2_n1 = 1, two.tailed = TRUE, var.equal = TRUE,
					 means=NULL, sds=NULL, conf.level = .95,
					 gen_fun=gen_t.test, return_analysis = FALSE, ...) {
	if(is.null(means))
		if(!missing(d) && !is.null(r))
			stop('Please use either d or r')
	if(!is.null(r)){
		type <- 'two.sample'
		stopifnot(var.equal)
	}
	dat <- gen_fun(n=n, n2_n1=n2_n1, d=d, r=r, type=type,
				   means=means, sds=sds, ...)
	res <- if(type == 'paired'){
		t.test(dat[[1]], dat[[2]], mu=mu,
			   paired=TRUE, conf.level=conf.level)
	} else if(type == 'two.sample'){
		t.test(dat[[1]], dat[[2]], var.equal=var.equal,
			   mu=mu, conf.level=conf.level)
	} else if(type == 'one.sample') {
		t.test(dat[[1]], mu=mu, conf.level=conf.level)
	}
	if(return_analysis) return(res)
	p <- res$p.value
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' @rdname p_t.test
#' @export
gen_t.test <- function(n, d, n2_n1 = 1, r = NULL, type = 'two.sample',
					   means=NULL, sds=NULL, ...){
	if(!is.null(r)){
		type <- 'two.sample'
		d <- r2d(r, n0=n, n1=n*n2_n1)
	}
	if(type == 'paired'){
		# n <- n * 2
		n2_n1 <- 1
	}
	n.each <- n * n2_n1
	# stopifnot(all.equal(n.each, as.integer(n.each)))
	if(type == 'one.sample'){
		DV <- if(!is.null(means))
			rnorm(n, mean=means, sd=sds) else rnorm(n, mean=d)
		dat <- list(DV)
	} else {
		if(!is.null(means)){
			if(!missing(d)) stop('d argument cannot be used with raw_info')
			group1 <- rnorm(n, mean=means[1], sd=sds[1])
			group2 <- rnorm(n * n2_n1, mean=means[2], sd=sds[2])
		} else {
			if(is.null(sds)) sds <- c(1,1)
			group1 <- rnorm(n, mean=d, sd = sds[1])
			group2 <- rnorm(n * n2_n1, sd = sds[2])
		}
		dat <- list(group1, group2)
	}
	dat
}

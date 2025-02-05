#' p-value from variance test simulation
#'
#' Generates one or or more sets of continuous data group-level data
#' to perform a variance test, and return a p-value. When two-samples
#' are investigated the \code{\link{var.test}} function will be used,
#' otherwise functions from the \code{EnvStats} package will be used.
#'
#'
#' @param n sample size per group, assumed equal across groups
#' @param sds a vector of standard deviations to use for each group
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param test type of test to use in multi-sample applications.
#'   Can be either \code{'Levene'} (default) or \code{'Bartlett'}
#' @param sigma standard deviation value to test against in one-sample test
#' @param correct logical; use correction when \code{test = 'Bartlett'}?
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{matrix} with k rows and k columns
#'   of counts. Default uses \code{\link{gen_var.test}}.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#'
#' @seealso \code{\link{gen_var.test}}
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
					   sigma = 1, two.tailed = TRUE,
					   test = 'Levene', correct = TRUE,
					   gen_fun=gen_var.test, ...){
	dat <- gen_fun(n=n, n.ratios=n.ratios, sds=sds, ...)
	p <- if(length(sds) == 1){
		EnvStats::varTest(dat$DV, sigma.squared = sigma^2)$p.value
	} else {
		if(length(sds) == 2) with(dat, var.test(DV ~ group)$p.value)
		else with(dat,
				  EnvStats::varGroupTest(DV ~ group, test=test, correct=correct)$p.value)
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' @rdname p_var.test
#' @export
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

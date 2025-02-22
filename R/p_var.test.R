#' p-value from variance test simulation
#'
#' Generates one or or more sets of continuous data group-level data
#' to perform a variance test, and return a p-value. When two-samples
#' are investigated the \code{\link{var.test}} function will be used,
#' otherwise functions from the \code{EnvStats} package will be used.
#'
#'
#' @param n sample size per group, assumed equal across groups
#' @param vars a vector of variances to use for each group; length of 1 for
#'   one-sample tests
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param test type of test to use in multi-sample applications.
#'   Can be either \code{'Levene'} (default) or \code{'Bartlett'}
#' @param sigma2 population variance to test against in one-sample test
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
#' p_var.test(100, vars=10, sigma2=9)
#'
#' # three sample
#' p_var.test(100, vars=c(10, 9, 11))
#'
#' if(FALSE){
#'   # power to detect three-group variance differences
#'   p_var.test(n=100, vars=c(10,9,11)) |> Spower()
#'
#'   # sample size per group to achieve 80% power
#'   p_var.test(n=NA, vars=c(10,9,11)) |>
#'          Spower(power=.80, interval=c(100, 1000))
#' }
#'
p_var.test <- function(n, vars, n.ratios = rep(1, length(vars)),
					   sigma2 = 1, two.tailed = TRUE,
					   test = 'Levene', correct = TRUE,
					   gen_fun=gen_var.test, ...){
	dat <- gen_fun(n=n, n.ratios=n.ratios, vars=vars, ...)
	p <- if(length(vars) == 1){
		EnvStats::varTest(dat$DV, sigma.squared = sigma2)$p.value
	} else {
		if(length(vars) == 2) with(dat, var.test(DV ~ group)$p.value)
		else with(dat,
				  EnvStats::varGroupTest(DV ~ group, test=test, correct=correct)$p.value)
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' @rdname p_var.test
#' @export
gen_var.test <- function(n, vars, n.ratios = rep(1, length(vars)), ...){
	sds <- sqrt(vars)
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

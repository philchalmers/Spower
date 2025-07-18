.SpowerEnv <- new.env(parent = emptyenv())

#' Get previously evaluated Spower execution
#'
#' If the result of \code{\link{Spower}} or \code{\link{SpowerBatch}} was not stored into
#' an object this function will retrieve the last evaluation.
#'
#' @return the last object returned from \code{\link{Spower}} or \code{\link{SpowerBatch}}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export
#'
getLastSpower <- function() .SpowerEnv$lastSim

#' Evaluate whether parameter is outside a given confidence interval
#'
#' Returns \code{TRUE} if parameter reflecting a null hypothesis
#' falls outside a given confidence interval. This is an alternative approach
#' to writing an experiment that returns a p-value.
#'
#' @param P0 parameter to evaluate
#' @param CI confidence interval
#' @return logical
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{is.CI_within}}, \code{\link{Spower}}
#' @export
#' @examples
#'
#' p0 <- .3
#' CI <- c(.2, .4)
#' is.outside_CI(p0, CI)
#'
#' # complement indicates if p0 is within CI
#' !is.outside_CI(p0, CI)
#'
#'
is.outside_CI <- function(P0, CI){
	stopifnot(length(P0) == 1 || length(CI) == 2)
	P0 < CI[1] || CI[2] < P0
}

#' Evaluate whether a confidence interval is within a tolerable interval
#'
#' Return \code{TRUE} if an estimated confidence interval falls within
#' a tolerable \code{interval} range. Typically used for
#' equivalence, superiority, or non-inferiority testing.
#'
#' @param CI estimated confidence interval (length 2)
#' @param interval tolerable interval range (length 2)
#' @return logical
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{is.outside_CI}}, \code{\link{Spower}}
#' @export
#'
#' @examples
#'
#' CI <- c(.2, .4)
#' LU <- c(.1, .3)
#' is.CI_within(CI, LU)        # not within tolerable interval
#' is.CI_within(CI, c(0, .5))  # is within wider interval
#'
#' # complement indicates if CI is outside interval
#' !is.CI_within(CI, LU)
#'
#' #####
#' # for superiority test
#' is.CI_within(CI, c(.1, Inf))  # CI is within tolerable interval
#'
#' # for inferiority test
#' is.CI_within(CI, c(-Inf, .3))  # CI is not within tolerable interval
#'
is.CI_within <- function(CI, interval){
	stopifnot(length(interval) == 2 || length(CI) == 2)
	interval[1] < CI[1] && CI[2] < interval[2]
}

Internal_Summarise <- function(condition, results, fixed_objects) {
	if(!is.null(fixed_objects$select))
		results <- results[ ,fixed_objects$select, drop=FALSE]
	results <- as.matrix(results)
	ret <- c(power = SimDesign::EDR(results[,1], alpha = condition$sig.level, unname=TRUE))
	ret
}

Internal_Summarise.Full <- function(condition, results, fixed_objects) {
	if(!is.null(fixed_objects$select))
		results <- results[ ,fixed_objects$select, drop=FALSE]
	results <- as.matrix(results)
	ret<- c(power = SimDesign::EDR(results, alpha = condition$sig.level,
						unname=ifelse(ncol(results) > 1, FALSE, TRUE)))
	ret
}

Internal_Summarise4Compromise <- function(condition, results, fixed_objects = NULL) {
	if(!is.null(fixed_objects$select))
		results <- results[ ,fixed_objects$select, drop=FALSE]
	results <- as.matrix(results)
	rate <- SimDesign::EDR(results[,1], alpha=condition$sig.level)
	ret <- c(beta_alpha = unname((1-rate) / condition$sig.level))
	ret
}

has.decimals <- function(x){
	intx <- as.integer(x)
	isTRUE(any(abs(x - intx) > 0))
}

# compute beta/alpha ratio given different alpha
compromise <- function(alpha, sim, Design, Summarise){
	Design$sig.level <- alpha
	out <- SimDesign::reSummarise(Summarise, results=sim, Design=Design)
	out$beta_alpha
}

compromise_root <- function(alpha, beta_alpha, ...)
	compromise(alpha, ...) - beta_alpha


# convert r to d
r2d <- function(rho, n0, n1){
	N <- n0 + n1
	d <- (N*rho) / sqrt(n0*n1*(1-rho^2))
	d
}

parent_env_fun <- function(level=2){
	ret <- NULL
	for(lev in level:2){
		nms <- ls(envir = parent.frame(lev))
		ret <- c(ret, nms)
	}
	ret
}

clip_CI <- function(CI){
	CI[CI > 1] <- 1
	CI[CI < 0] <- 0
	CI
}

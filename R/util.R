.SpowerEnv <- new.env(parent = emptyenv())

#' Get previously evaluated Spower execution
#'
#' If the result of \code{\link{Spower}} was not stored into
#' an object this function will retrieve the last evaluation.
#'
#' @return the last object returned from \code{\link{Spower}}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{Spower}}
#' @export
#'
getLastSpower <- function() .SpowerEnv$lastSim

#' Evaluate whether parameter is outside a given confidence interval
#'
#' Return a logical if parameter reflecting a null hypothesis
#' falls outside a given confidence interval.
#'
#' @param P0 parameter to evaluate
#' @param CI confidence interval
#' @return logical
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{is.within_CI}}
#' @export
#' @examples
#'
#' p0 <- .3
#' CI <- c(.2, .4)
#' is.outside_CI(p0, CI)
#' is.within_CI(p0, CI)
#'
#'
is.outside_CI <- function(P0, CI){
	P0 < CI[1] || CI[2] < P0
}

#' Evaluate whether parameter is within a given confidence interval
#'
#' Return a logical if parameter reflecting a null hypothesis
#' falls within a given confidence interval.
#'
#' @param P0 parameter to evaluate
#' @param CI confidence interval
#' @return logical
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{is.outside_CI}}
#' @export
#'
#' @examples
#'
#' p0 <- .3
#' CI <- c(.2, .4)
#' is.outside_CI(p0, CI)
#' is.within_CI(p0, CI)
#'
is.within_CI <- function(P0, CI){
	CI[1] < P0 && P0 < CI[2]
}

Internal_Summarise <- function(condition, results, fixed_objects) {
	results <- as.matrix(results)
	ret <- c(power = EDR(results[,1], alpha = condition$sig.level, unname=TRUE))
	ret
}

Internal_Summarise.Full <- function(condition, results, fixed_objects) {
	results <- as.matrix(results)
	ret<- c(power = EDR(results, alpha = condition$sig.level,
						unname=ifelse(ncol(results) > 1, FALSE, TRUE)))
	ret
}

Internal_Summarise4Compromise <- function(condition, results, fixed_objects = NULL) {
	results <- as.matrix(results)
	rate <- EDR(results[,1], alpha=condition$sig.level)
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
	out <- reSummarise(Summarise, results=sim, Design=Design)
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

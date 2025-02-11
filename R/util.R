.SpowerEnv <- new.env(parent = emptyenv())

#' Get previously evaluated Spower execution
#'
#' If the result of \code{\link{Spower}} was not stored into
#' an object this function will retrieve the last evaluation.
#'
#' @return the last object returned from \code{\link{Spower}}
#' @seealso \code{\link{Spower}}
#' @export
#'
getLastSpower <- function() .SpowerEnv$lastSim

Internal_Summarise <- function(condition, results, fixed_objects) {
	ret <- if(is.null(dim(results)))
		c(power = EDR(results, alpha = condition$sig.level))
	else {
		c(power = EDR(results[,1], alpha = condition$sig.level))
	}
	ret
}

Internal_Summarise4Compromise <- function(condition, results, fixed_objects = NULL) {
	rate <- EDR(results, alpha=condition$sig.level)
	ret <- c(beta_alpha = unname((1-rate) / condition$sig.level))
	ret
}

Internal_Summarise4Compromise.empirical <- function(condition,
													results, fixed_objects = NULL) {
	rates <- EDR(results, alpha=condition$sig.level)
	ret <- c(beta_alpha = unname((1-rates['power']) / rates['alpha']))
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
		is_fun <- sapply(nms, function(x, envir) is.function(get(x, envir=envir)),
						 envir = parent.frame(lev))
		if(any(is_fun)) ret <- c(ret, nms[is_fun])
	}
	ret
}

Internal_Summarise <- function(condition, results, fixed_objects) {
	ret <- c(power = EDR(results, alpha = condition$sig.level))
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

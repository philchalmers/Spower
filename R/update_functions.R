#' Update compromise analysis criteria without re-simulating
#'
#' When a power or compromise analysis was performed in
#' \code{\link{Spower}} this function can be used to
#' update the compromise criteria without the need for re-simulating
#' the experiment.
#'
#' @param x object returned from \code{\link{Spower}} where \code{power}
#'   was estimated or the \code{bete_alpha} criteria were supplied
#' @param beta_alpha Type II/Type I error ratio
#' @export
#'
#' @examples
#' \dontrun{
#' # Solve beta/alpha ratio to specific error trade-off constant
#' (out <- Spower(p_t.test, n = 50, d = .5, beta_alpha = 2))
#'
#' # update beta_alpha criteria without re-simulating
#' updateCompromise(out, beta_alpha=4)
#'
#' }
#'
updateCompromise <- function(x, beta_alpha){
	ret <- x
	conditions <- attr(x, 'Spower_extra')$conditions
	conditions$beta_alpha <- NULL
	out <- uniroot(compromise_root, c(.0001, .9999), beta_alpha=beta_alpha,
				   sim=ret, Design=conditions, Summarise=Internal_Summarise4Compromise)
	ret$sig.level <- out$root
	ret$power <- 1 - beta_alpha * out$root
	conditions$sig.level <- as.numeric(NA)
	conditions$beta_alpha <- beta_alpha
	attr(ret, 'Spower_extra')$conditions <- conditions
	ret
}

#' Update power estimates given new sig.level without re-simulating
#'
#' When a power was performed with \code{\link{Spower}} this
#' function can be used to
#' update the power estimate given an alternative Type I error rate
#' (alpha/sig.level) without re-simulating the experiment.
#'
#' @param x object returned from \code{\link{Spower}} where \code{power}
#'   was estimated or the \code{bete_alpha} criteria were supplied
#' @param sig.level Type I error rate (alpha)
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # estimate power using sig.level = .05 (default)
#' (out <- Spower(p_t.test, n = 50, d = .5))
#'
#' # update power estimate given sig.level = .01
#' update_sig.level(out, sig.level=.01)
#'
#' }
#'
update_sig.level <- function(x, sig.level){
	if(missing(sig.level))
		stop('Must specify sig.level')
	stopifnot(sig.level < 0 || sig.level > 1)
	browser()
}

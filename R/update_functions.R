#' Update compromise or post-hoc power analysis without re-simulating
#'
#' When a power or compromise analysis was performed in
#' \code{\link{Spower}} this function can be used to
#' update the compromise or power criteria without the need for re-simulating
#' the experiment. For compromise analyses a \code{beta_alpha} criteria
#' must be supplied, while for post-hoc power analyses the \code{sig.level}
#' must be supplied.
#'
#' @param object object returned from \code{\link{Spower}} where \code{power}
#'   was estimated or the \code{bete_alpha} criteria were supplied
#' @param beta_alpha Type II/Type I error ratio
#' @param sig.level Type I error rate (alpha)
#' @param predCI confidence interval precision (see \code{\link{Spower}} for
#'   similar input)
#' @param ... arguments to be passed
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ########
#' ## Post-hoc power analysis update
#'
#' # Estimate power using sig.level = .05 (default)
#' out <- Spower(p_t.test, n = 50, d = .5)
#'
#' # update power estimate given sig.level=.01 and .20
#' update(out, sig.level=.01)
#' update(out, sig.level=.20)
#'
#'
#' ########
#' ## Compromise analysis update
#'
#' # Solve beta/alpha ratio to specific error trade-off constant
#' out <- Spower(p_t.test, n = 50, d = .5, beta_alpha = 2)
#'
#' # update beta_alpha criteria without re-simulating
#' update(out, beta_alpha=4)
#'
#' }
#'
#' @importFrom stats update


#' @rdname update
#' @export
update.Spower <- function(object, sig.level = .05, beta_alpha = NULL, predCI=.95, ...){
	ret <- if(!is.null(beta_alpha)){
		updateCompromise(object, beta_alpha=beta_alpha)
	} else {
		update_sig.level(object, sig.level=sig.level, predCI=predCI)
	}
	ret
}

updateCompromise <- function(x,  beta_alpha = NULL){
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

update_sig.level <- function(x, sig.level, predCI=.95){
	if(missing(sig.level))
		stop('Must specify sig.level')
	stopifnot(sig.level > 0 && sig.level < 1)
	x$sig.level <- attr(x, "Spower_extra")$conditions$sig.level <- sig.level
	out <- SimDesign::reSummarise(Internal_Summarise, results=x)
	design_names <- attr(x, 'design_names')$design
	x$power <- as.numeric(out[,length(design_names)+1])
	replications <- x$REPLICATIONS
	alpha <- 1 - predCI
	CI <- x$power + c(qnorm(c(alpha/2, predCI+alpha/2))) *
		sqrt((x$power * (1-x$power))/replications)
	CI[CI > 1] <- 1
	CI[CI < 0] <- 0
	names(CI) <- paste0('CI_', c(alpha/2, predCI+alpha/2)*100)
	attr(x, 'extra_info')$power.CI <- CI
	x
}

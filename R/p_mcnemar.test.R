#' p-value from McNemar test simulation
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
#' @param OR instead of supplying the \code{prop} table, the
#'   odds ratio can be specified instead \eqn{\pi_{12}/\pi_{21}}. Also requires
#'   proportion of discordant pairings to be specified
#' @param prop.disc proportion of discordant pairings, \eqn{\pi_{12} + \pi_{21}}
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{matrix} with k rows and k columns
#'   of counts. Default uses \code{\link{gen_mcnemar.test}}.
#'   User defined version of this function must include the argument \code{...}
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#'
#' @seealso \code{\link{gen_mcnemar.test}}
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
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
#' # return analysis model
#' p_mcnemar.test(n=sum(Performance), prop=prop, return_analysis=TRUE)
#'
#' \donttest{
#'
#' # post-hoc power (not recommended)
#' Spower(p_mcnemar.test(n=sum(Performance), prop=prop))
#'
#' # odds ratio + discordant proportions supplied instead
#' OR <- prop[1,2] / prop[2,1]
#' disc <- prop[1,2] + prop[2,1]
#' p_mcnemar.test(n=50, OR=.25, prop.disc=disc, two.tailed=FALSE) |>
#'   Spower(replications=30000)
#'
#' }
#'
#' @export
p_mcnemar.test <- function(n, prop, OR=NULL, prop.disc=NULL,
						   two.tailed = TRUE, correct=TRUE,
						   gen_fun=gen_mcnemar.test, return_analysis = FALSE, ...) {
	if(!is.null(OR)){
		prop <- matrix(0, 2, 2)
		prop[2,1] <- prop.disc / (OR + 1)
		prop[1,2] <- prop.disc - prop[2,1]
		prop[1,1] <- prop[2,2] <- (1 - sum(prop))/2
	}
	dat <- gen_fun(n=n, prop=prop, ...)
	ret <- mcnemar.test(dat, correct=correct)
	if(return_analysis) return(ret)
	p <- ret$p.value
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' @rdname p_mcnemar.test
#' @export
gen_mcnemar.test <- function(n, prop, ...) {
	draws <- rmultinom(1, n, prob = as.numeric(prop))
	dat <- matrix(draws, nrow(prop), ncol(prop))
	dat
}

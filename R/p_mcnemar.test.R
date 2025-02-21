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
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{matrix} with k rows and k columns
#'   of counts. Default uses \code{\link{gen_mcnemar.test}}.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#'
#' @seealso \code{\link{gen_mcnemar.test}}
#' @return a single p-value
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
#' \dontrun{
#'
#' # post-hoc power (not recommended)
#' Spower(p_mcnemar.test(n=sum(Performance), prop=prop))
#'
#' }
#'
#' @export
p_mcnemar.test <- function(n, prop,
						   two.tailed = TRUE, correct=TRUE,
						   gen_fun=gen_mcnemar.test, ...) {
	dat <- gen_fun(n=n, prop=prop, ...)
	p <- mcnemar.test(dat, correct=correct)$p.value
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

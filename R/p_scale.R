#' p-value from Scale Test simulation
#'
#' Simulates data given one or two parent distributions and
#' returns a p-value testing that the scale of the type distributions are the same.
#' Default implementation uses Gaussian distributions, however the
#' distribution function may be modified to
#' reflect other populations of interest.
#' Uses \code{\link{ansari.test}} or \code{\link{mood.test}} for the analysis.
#'
#' @param n sample size per group
#' @param scale the scale to multiply the second group by (1 reflects equal scaling)
#' @param n2_n1 sample size ratio
#' @param test type of method to use. Can be either \code{'Ansari'} or \code{'Mood'}
#' @param exact a logical indicating whether an exact p-value should be computed
#' @param parent data generation function (default assumes Gaussian shape). Must be
#'   population mean centered
#' @param two.tailed logical; use two-tailed test?
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @param ... additional arguments to pass to simulation functions (if used)
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export
#' @examples
#'
#' # n=30 per group,
#' #  Distributions Gaussian with sd=1 for first group and sd=2 for second
#' p_scale(30, scale=2)
#' p_scale(30, scale=2, test='Mood')
#'
#' # compare chi-squared distributions
#' parent <- function(n, df, ...) rchisq(n, df=df) - df
#' p_scale(30, scale=2, parent=parent, df=3)
#'
#' \donttest{
#'   # empirical power of the experiments
#'   p_scale(30, scale=2) |> Spower()
#'   p_scale(30, scale=2, test='Mood') |> Spower()
#'
#'   p_scale(30, scale=2, parent=parent, df=3) |> Spower()
#'   p_scale(30, scale=2, test='Mood', parent=parent, df=3) |> Spower()
#'
#' }
p_scale <- function(n, scale, n2_n1 = 1, two.tailed = TRUE,
					exact = NULL, test = 'Ansari',
					parent = function(n, ...) rnorm(n), ...,
					return_analysis = FALSE){
	stopifnot(test %in% c('Ansari', 'Mood'))
	dat1 <- parent(n, ...)
	dat2 <- parent(n*n2_n1, ...) * scale
	ret <- if(test == 'Ansari')
		ansari.test(dat1, dat2, exact=exact, alternative='two.sided')
	else mood.test(dat1, dat2, alternative='two.sided')
	if(return_analysis) return(ret)
	p <- ret$p.value
	p <- ifelse(two.tailed, p, p/2)
	p
}

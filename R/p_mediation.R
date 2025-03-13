#' p-value from three-variable mediation analysis simulation
#'
#' Simple 3-variable mediation analysis simulation to test the hypothesis that
#' X -> Y is mediated by the relationship X -> M -> Y. Currently,
#' M and Y are assumed to be continuous variables with Gaussian errors, while
#' X may be continuous or dichotomous.
#'
#' @param n total sample size unless \code{dichotomous.X = TRUE}, in which the
#'   value represents the size per group
#' @param a regression coefficient for the path X -> M
#' @param b regression coefficient for the path M -> Y
#' @param cprime partial regression coefficient for the path X -> Y
#' @param dichotomous.X logical; should the X variable be generated as though it
#'  were dichotomous? If TRUE then \code{n} represents the sample size per group
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param method type of inferential method to use. Default uses the Wald
#'    (a.k.a., Sobel) test
#' @param gen_fun function used to generate the required two-sample data.
#'   Object returned must be a \code{data.frame} with the columns
#'   \code{"DV"} and \code{"group"}. Default uses \code{\link{gen_mediation}}
#'   to generate conditionally Gaussian distributed samples.
#'   User defined version of this function must include the argument \code{...}
#' @param sd.X standard deviation for X
#' @param sd.M standard deviation for M
#' @param sd.Y standard deviation for Y
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#'
#' @seealso \code{\link{gen_mediation}}
#' @importFrom lavaan sem lavInspect parameterEstimates
#' @return a single p-value
#' @examples
#'
#' # joint test H0: a*b = 0
#' p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39)
#' p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, dichotomous.X=TRUE)
#'
#' \dontrun{
#'
#'   # power to detect mediation
#'   p_mediation(n=50, a=sqrt(.35), b=sqrt(.35), cprime=.39) |>
#'     Spower(parallel=TRUE, replications=1000)
#'
#'   # sample size estimate for .95 power
#'   p_mediation(n=NA, a=sqrt(.35), b=sqrt(.35), cprime=.39) |>
#'     Spower(power=.95, interval=c(50, 200), parallel=TRUE)
#'
#' }
#'
#' @export
p_mediation <- function(n, a, b, cprime, dichotomous.X=FALSE,
						two.tailed=TRUE, method = 'wald',
						sd.X=1, sd.Y=1, sd.M=1,
						gen_fun=gen_mediation, ...){
	dat <- gen_fun(n, a=a, b=b, cprime=cprime,
				   sd.X=sd.X, sd.Y=sd.Y, sd.M=sd.M,
				   dichotomous.X=dichotomous.X, ...)
	model <- '
	       # direct effect
             Y ~ c*X
           # mediator
             M ~ a*X
             Y ~ b*M
           # indirect effect (a*b)
             ab := a*b
         '
	fit <- lavaan::sem(model, data=dat)
	if(!lavaan::lavInspect(fit, 'converged'))
		stop('Model did not converge')
	PE <- lavaan::parameterEstimates(fit)
	if(method == 'wald')
		p <- PE$pvalue[PE$lhs == 'ab']   # joint test
	else stop('other methods not yet supported')
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' @rdname p_mediation
#' @export
gen_mediation <- function(n, a, b, cprime, dichotomous.X=FALSE,
						  sd.X=1, sd.Y=1, sd.M=1, ...){
	if(dichotomous.X){
		X <- rep(0:1, each=n)
		X <- X * sd.X / sqrt(.25)
		n <- n*2
	} else X <- rnorm(n)
	M <- a*X + rnorm(n, sd=sd.M - a^2)
	Y <- b*M + cprime*X + rnorm(n, sd=sd.Y - b^2 - cprime^2)
	dat <- data.frame(X=X * sd.X, Y=Y * sd.Y, M=M * sd.M)
	dat
}

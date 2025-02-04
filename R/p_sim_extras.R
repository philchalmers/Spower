#' p-value from three-variable mediation analysis simulation
#'
#' Simple 3-variable mediation analysis simulation to test the hypothesis that
#' X -> Y is mediated by the relationship X -> M -> Y. Currently,
#' M and Y are assumed to be continuous variables with Gaussian errors, while
#' X may be continuous or dichotomous.
#'
#' @param n total sample size
#' @param a regression coefficient for the path X -> M
#' @param b regression coefficient for the path M -> Y
#' @param cprime partial regression coefficient for the path X -> Y
#' @param dichotomous.X logical; should the X variable be generated as though it
#'  were dichotomous? If TRUE then \code{n} represents the sample size per group
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param method type of inferential method to use. Default uses the Wald
#'    (a.k.a., Sobel) test
#' @importFrom lavaan sem lavInspect parameterEstimates
#' @return a single p-value
#' @examples
#'
#' # joint test H0: a*b = 0
#' p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39)
#' p_mediation(50, a=sqrt(.35), b=sqrt(.35), cprime=.39, dichotomous.X=TRUE)
#'
#' if(FALSE){
#'
#'   # power to detect mediation
#'   Spower(p_mediation, n=50, a=sqrt(.35), b=sqrt(.35),
#'          cprime=.39, parallel=TRUE, replications=1000)
#'
#'   # sample size for .95 power
#'   Spower(p_mediation, n=NA, a=sqrt(.35), b=sqrt(.35),
#'          cprime=.39, power=.95, interval=c(50, 300),
#'          parallel=TRUE)
#'
#'
#' }
#'
#' @export
p_mediation <- function(n, a, b, cprime, dichotomous.X=FALSE,
						two.tailed=TRUE, method = 'wald'){
	if(dichotomous.X){
		X <- rep(0:1, each=n)
		n <- n*2
	} else X <- rnorm(n)
	M <- a*X + rnorm(n)
	Y <- b*M + cprime*X + rnorm(n)
	dat <- data.frame(X, Y, M)
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




#' p-value from tetrachoric/polychoric or polyserial
#'
#' Generates correlated X-Y data and returns a p-value to assess the null
#' of no correlation in the population. The X-Y data are generated
#' assuming a multivariate normal distribution and subsequently discretized
#' for one or both of the variables.
#'
#' @param n sample size
#' @param r correlation prior to the discretization (recovered via the
#'   polyserial/polychoric estimates)
#' @param rho population coefficient to test against
#' @param tauX intercept parameters used for discretizing the X variable
#' @param tauY intercept parameters used for discretizing the Y variable. If
#'   missing a polyserial correlation will be estimated, otherwise a
#'   tetrachoric/polychoric correlation will be estimated
#' @param ML logical; use maximum-likelihood estimation?
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param score logical; should the SE be based at the null hypothesis (score test)
#'   or the ML estimate (Wald test)? The former is the canonical form for
#'   a priori power analyses though requires twice as many computations as the
#'   Wald test approach
#' @param gen_fun function used to generate the required
#'   continuous bivariate data (prior to truncation).
#'   Object returned must be a \code{matrix} with two columns.
#'   Default uses \code{\link{gen_r}} to generate conditionally
#'   dependent data from a bivariate normal distribution.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#'
#' @seealso \code{\link{gen_r}}
#' @return a single p-value
#' @export
#' @examples
#'
#' # 100 observations, .5 correlation, tetrachoric estimate
#' p_r.cat(100, r=.5, tauX=0, tauY=1)
#'
#' # Wald test
#' p_r.cat(100, r=.5, tauX=0, tauY=1, score=FALSE)
#'
#' # polyserial estimate (Y continuous)
#' p_r.cat(50, r=.5, tauX=0)
#'
p_r.cat <- function(n, r, tauX, rho=0, tauY = NULL,
					ML=TRUE, two.tailed=TRUE, score=FALSE,
					gen_fun=gen_r, ...){
	continuous.Y <- is.null(tauY)
	dat <- gen_fun(n=n, r=r, ...)
	datcut <- matrix(0, n, 2)
	for(i in length(tauX):1)
		datcut[dat[,1] > tauX[i], 1] <- i
	if(!continuous.Y){
		for(i in length(tauY):1)
			datcut[dat[,2] > tauY[i], 2] <- i
	} else datcut[,2] <- dat[,2]
	datcut <- as.data.frame(datcut)
	colnames(datcut) <- c('x', 'y')
	# could generate r=0 data to get SE_0 instead for proper score test
	out <- if(continuous.Y){
		with(datcut, polycor::polyserial(y, x, ML=ML, std.err=TRUE))
	} else {
		with(datcut, polycor::polychor(y, x, ML=ML, std.err=TRUE))
	}
	est <- out$rho
	vcov <- out$var
	if(score > 1) return(sqrt(vcov[1,1]))
	SE <- if(score == 1){
		p_r.cat(n, r=rho, tauX=tauX, rho=rho, tauY=tauY,
				ML=ML, two.tailed=two.tailed, score=2)
	} else sqrt(vcov[1,1])
	z <- (est - rho) / SE
	p <- pnorm(abs(z), lower.tail=FALSE)
	p <- ifelse(two.tailed, p*2, p)
	p
}

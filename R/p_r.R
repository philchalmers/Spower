#' p-value from correlation simulation
#'
#' Generates correlated X-Y data and returns a p-value to assess the null
#' of no correlation in the population. The X-Y data are generated
#' assuming a bivariate normal distribution.
#'
#' @param n sample size
#' @param r correlation
#' @param rho population coefficient to test against. Uses the
#'   Fisher's z-transformation approximation when non-zero
#' @param method method to use to compute the correlation
#'   (see \code{\link{cor.test}}). Only used when \code{rho = 0}
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param gen_fun function used to generate the required dependent bivariate data.
#'   Object returned must be a \code{matrix} with two columns and \code{n} rows.
#'   Default uses \code{\link{gen_r}} to generate conditionally
#'   dependent data from a bivariate normal distribution.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#' @seealso \code{\link{gen_r}}
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @examples
#'
#' # 50 observations, .5 correlation
#' p_r(50, r=.5)
#' p_r(50, r=.5, method = 'spearman')
#'
#' # test against constant other than rho = .6
#' p_r(50, .5, rho=.60)
#'
#'
#' \donttest{
#'     # compare simulated results to pwr package
#'
#'     pwr::pwr.r.test(r=0.3, n=50)
#'     p_r(n=50, r=0.3) |> Spower()
#'
#'     pwr::pwr.r.test(r=0.3, power=0.80)
#'     p_r(n=NA, r=0.3, power=.80, interval=c(10, 200)) |> Spower()
#'
#'     pwr::pwr.r.test(r=0.1, power=0.80)
#'     p_r(n=NA, r=0.1, power=.80, interval=c(200, 1000)) |> Spower()
#'
#' }
#'
#' @export
p_r <- function(n, r, rho = 0, method = 'pearson', two.tailed = TRUE,
				gen_fun=gen_r, ...) {
	dat <- gen_fun(n=n, r=r, ...)
	colnames(dat) <- c('x', 'y')
	if(rho != 0){
		out <- cor.test(~ x + y, dat, method=method)
		z <- with(out, 1/2 * log((1+estimate)/(1-estimate)))
		se <- 1 / sqrt(n-3)
		z0 <- 1/2 * log((1+rho)/(1-rho))
		t <- (z - z0) / se
		p <- pnorm(abs(t), lower.tail=FALSE)*2
	} else {
		p <- cor.test(~ x + y, dat, method=method)$p.value
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' @rdname p_r
#' @export
gen_r <- function(n, r, ...){
	dat <- SimDesign::rmvnorm(n, sigma = matrix(c(1,r,r,1), 2, 2))
	dat
}

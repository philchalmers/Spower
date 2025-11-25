#' p-value from simple linear regression model simulation
#'
#' p-values associated with the simple linear regression model,
#' \eqn{y = \beta_0 + \beta_1 X + \epsilon}.
#' Focus is on the slope/intercept behavior of the model.
#'
#' @param n sample size
#' @param beta slope parameter
#' @param beta0 null value to test against
#' @param test test to evaluate using \code{\link[car]{lht}}.
#'   Default evaluates the null hypothesis that the slope is equal to 0
#' @param sd_x standard deviation of IV
#' @param sd_y standard deviation of DV
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @param gen_fun function used to generate the required X-Y data.
#'   Object returned must be a \code{data.frame} with the columns
#'   \code{'y'} and \code{'x'}. Default uses \code{\link{gen_slr}}.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#' @seealso \code{\link{p_glm}}, \code{\link{p_lm.R2}}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @return a single p-value
#' @export
#' @examples
#'
#' p_slr(n=100, beta = -0.0667, sd_x = 7.5, sd_y = 4)
#'
#' \dontrun{
#' p_slr(n=100, beta = -0.0667, sd_x = 7.5, sd_y = 4) |> Spower()
#' }
#'
p_slr <- function(n, beta, sd_x, sd_y, beta0 = 0, test = 'x = 0',
				  gen_fun=gen_slr,
					return_analysis = FALSE, ...){
	df <- gen_fun(n=n, beta=beta, sd_x=sd_x, sd_y=sd_y, ...)
	mod <- lm(y ~ x, df)
	p <- car::lht(mod, test)$`Pr`[2]
	p
}

#' @rdname p_slr
#' @export
gen_slr <- function(n, beta, sd_x, sd_y, ...){
	rho <- beta * (sd_x/sd_y)
	rho_c <- rho * sd_x * sd_y
	dat <- SimDesign::rmvnorm(n=n,
			sigma=matrix(c(sd_y^2,rho_c, rho_c, sd_x^2), 2, 2))
	colnames(dat) <- c('y', 'x')
	dat <- as.data.frame(dat)
	dat
}

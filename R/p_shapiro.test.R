#' p-value from Shapiro-Wilk Normality Test simulation
#'
#' Generates univariate distributional data and returns a p-value to assess the null
#' that the population follows a Gaussian distribution shape. Uses
#' \code{\link{shapiro.test}}.
#'
#' @param dist expression used to generate the required sample data
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @examples
#'
#' # 50 observations drawn from normal distribution (null is true)
#' p_shapiro.test(rnorm(50))
#'
#' # return analysis object
#' p_shapiro.test(rnorm(50), TRUE)
#'
#' # 50 observations from slightly skewed chi-squared distribution (power)
#' p_shapiro.test(rchisq(50, df=100))
#'
#' \donttest{
#'     # empirical Type I error rate estimate
#'     p_shapiro.test(rnorm(50)) |> Spower()
#'
#'     # power
#'     p_shapiro.test(rchisq(50, df=100)) |> Spower()
#' }
#'
#' @export
p_shapiro.test <- function(dist, return_analysis = FALSE) {
	res <- shapiro.test(dist)
	if(return_analysis) return(res)
	p <- res$p.value
	p
}


#' p-value from Shapiro-Wilk Normality Test simulation
#'
#' Generates univariate distributional data and returns a p-value to assess the null
#' that the population follows a Gaussian distribution shape. Uses
#' \code{\link{shapiro.test}}.
#'
#' @param n sample size
#' @param dist expression used to generate the required sample data
#' @return a single p-value
#' @examples
#'
#' # 50 observations drawn from normal distribution (null is true)
#' p_shapiro.test(rnorm(50))
#'
#' # 50 observations from slightly skewed chi-squared distribution (power)
#' p_shapiro.test(rchisq(50, df=100))
#'
#' \dontrun{
#'     # empirical Type I error rate estimate
#'     p_shapiro.test(rnorm(50)) |> Spower()
#'
#'     # power
#'     p_shapiro.test(rchisq(50, df=100)) |> Spower()
#' }
#'
#' @export
p_shapiro.test <- function(dist) {
	p <- shapiro.test(dist)$p.value
	p
}


#' p-value from Mauchly's Test of Sphericity simulation
#'
#' Perform simulation experiment for Mauchly's Test of Sphericity using
#' the function \code{mauchlys.test}, returning a p-value.
#' Assumes the data are from a multivariate
#' normal distribution, however this can be modified.
#'
#' @param n sample size
#' @param sigma symmetric covariance/correlation matrix passed to \code{gen_fun}
#' @param gen_fun function used to generate the required data.
#'   Object returned must be a \code{matrix} with \code{K} columns and
#'   \code{n} rows. Default uses \code{\link{gen_mauchly.test}}
#'   to generate multivariate normal samples.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#' @return a single p-value
#' @examples
#'
#' sigma <- diag(c(1,2,1))
#' sigma
#'
#' p_mauchly.test(100, sigma=sigma)
#'
#' # Null is true
#' sigma.H0 <- diag(3)
#' p_mauchly.test(100, sigma=sigma.H0)
#'
#'
#' \dontrun{
#'     # empirical power estimate
#'     p_mauchly.test(100, sigma=sigma) |> Spower()
#'
#'     # empirical Type I error estimate
#'     p_mauchly.test(100, sigma=sigma.H0) |> Spower()
#' }
#'
#' @export
p_mauchly.test <- function(n, sigma, gen_fun=gen_mauchly.test, ...){
	dat <- gen_fun(n, sigma=sigma, ...)
	p <- mauchlys.test(dat)
	p
}

#' @rdname p_mauchly.test
#' @export
gen_mauchly.test <- function(n, sigma, ...){
	rmvnorm(n, sigma=sigma)
}

#' @rdname p_mauchly.test
#' @param X a matrix with \code{k} columns and \code{n} rows
#' @export
mauchlys.test <- function(X) {
	# code borrowed and modified from superb::MauchlySphericityTest,
	# version 0.95.9; date: 2025-03-13
	p   <- ncol(X)
	n   <- nrow(X)
	S   <- cov(X)
	Sij <- t(t(S - rowMeans(S)) - colMeans(S)) + mean(rowMeans(S))
	lam <- eigen(Sij)$values
	W   <- prod(lam[1:p-1])/(1/(p-1) * sum(lam[1:p-1]))^(p-1)
	f   <- (2 * (p - 1)^2 + p + 1)/(6 * (p - 1) * (n - 1))
	df  <- p * (p - 1)/2 - 1
	suppressWarnings(
		chiW <- if(is.na(log(W))) -(1 - f) * (n - 1) * log(abs(W))
		  else -(1 - f) * (n - 1) * log(W)
	)
	pW <- 1 - pchisq(chiW, df)
	pW
}

#' p-value from (generalized) linear regression model simulations with fixed predictors
#'
#' p-values associated with (generalized) linear regression model.
#' Requires a prespecified design matrix (\code{X}).
#'
#' @param formula formula passed to either \code{\link{lm}} or
#'   \code{\link{glm}}
#' @param X a data.frame containing the covariates
#' @param betas vector of slope coefficients that match the
#'   \code{model.matrix} version of \code{X}
#' @param test character vector specifying the test to pass to
#'   \code{\link[car]{lht}}. Can also be a list of character vectors
#'   to evaluate multiple tests
#' @param sigma residual standard deviation for linear model. Only
#'  used when \code{family = 'gaussian'}
#' @param family family of distributions to use (see \code{\link{family}})
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{data.frame}. Default uses \code{\link{gen_lm.R2}}.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#' @seealso \code{\link{gen_lm.R2}}
#' @export
#' @examples
#'
#' X <- data.frame(G = factor(rep(c('control', 'treatment'), each=50)),
#'                 C = sample(50:100, 100, replace=TRUE))
#' head(X)
#'
#' # ANCOVA setup
#' p_glm(y ~ G + C, test="Gtreatment = 0",
#'   X=X, betas=c(10, .3, 1), sigma=1)
#'
#' # ANCOVA setup with logistic regression
#' p_glm(y ~ G + C, test="Gtreatment = 0",
#'   X=X, betas=c(-2, .5, .01), family=binomial()) |> Spower()
#'
p_glm <- function(formula, X, betas, test, sigma = NULL,
				  family = gaussian(), gen_fun=gen_glm, ...){
	family <- add.sample2family(family)
	X <- gen_fun(formula=formula, X=X, betas=betas, sigma=sigma,
				 family=family, ...)
	mod <- if(family$family == 'gaussian'){
		lm(formula=formula, data=X, ...)
	} else {
		glm(formula=formula, data=X, family=family, ...)
	}
	if(!is.list(test))
		test <- list(test)
	ps <- sapply(test, \(tst) car::lht(mod, tst)$`Pr`[2])
	if(length(ps) > 1)
		names(ps) <- paste0('test_', 1:length(ps))
	ps
}

gen_glm <- function(formula, X, betas, sigma = NULL,
					family = gaussian(), ...){
	ynm <- as.character(formula[2])
	X[ynm] <- 0
	Xf <- model.frame(formula=formula, data=X)
	Xd <- model.matrix(formula, data=Xf)
	stopifnot(ncol(Xd) == length(betas))
	yhat <- betas %*% t(Xd)
	N <- length(yhat)
	if(family$family == 'gaussian'){
		stopifnot(!is.null(sigma))
		y <- yhat + rnorm(N, sd=sigma)
	} else {
		mus <- as.vector(family$linkinv(yhat))
		y <- family$sample(mus)
	}
	X[[ynm]] <- y
	X
}

add.sample2family <- function(family){
	if(family$family == 'gaussian') return(family)
	if(family$family == 'binomial'){
		family$sample <- function(mus)
			sapply(mus, \(mu) rbinom(1, size = 1, prob=mu))
	} else {
		stop('family data generator not yet supported')
	}
	family
}

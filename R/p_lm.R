#' p-value from global linear regression model simulation
#'
#' p-values associated with linear regression model using fixed or random
#' independent variables.
#'
#' @param n sample size
#' @param R2 R-squared effect size
#' @param k number of IVs
#' @param R2_0 null hypothesis for R-squared
#' @param k.R2_0 number of IVs associated with the null hypothesis model
#' @param R2.resid residual R-squared value, typically used when comparing
#'   nested models when fit sequentially (e.g., comparing model A vs B when
#'   model involves the structure A -> B -> C)
#' @param fixed.X logical; should the IVs be considered fixed or random?
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{data.frame}. Default uses \code{\link{gen_lm}}.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#' @seealso \code{\link{gen_lm}}
#' @export
#' @examples
#'
#' # 5 fixed IVs, R^2 = .1, sample size of 95
#' p_lm(n=95, R2=.1, k=5)
#'
p_lm <- function(n, R2, k, R2_0 = 0, k.R2_0 = 0,
				 R2.resid=1-R2, fixed.X=TRUE, gen_fun=gen_lm, ...){
	stopifnot(R2 > R2_0)
	if(k.R2_0 == 0)
		stopifnot(k >= 2)
	df <- gen_fun(n=n, fixed.X=fixed.X, k=k,
				  R2=R2, R2_0=R2_0, R2.resid=R2.resid, ...)
	if(k.R2_0 > 0)
		df2 <- df[,c(1,3:(k.R2_0 + 2))]
	mod1 <- lm(y ~ ., df)
	if(!fixed.X && k.R2_0 == 0 && R2_0 != 0){
		stop('Random X with non-zero R2_0 not currently supported', call.=FALSE)
	} else {
		mod0 <- if(R2_0 == 0) lm(y ~ 1, df)
		else lm(y ~ ., df2)
	}
	p <- anova(mod0, mod1)[2, "Pr(>F)"]
	p
}

if(FALSE){

	# Example 7.3b
	# G*power gives 0.3464 (broken)
	Spower(p_lm, n=100, R2=.4, R2_0 = .3, k=5, fixed.X=FALSE)

	# Example 7.3c
	# G*power gives N=153 (broken)
	Spower(p_lm, n=NA, R2=.05, R2_0 = .2, k=5, fixed.X=FALSE, power=.9,
		   interval=c(50,300))


}

#' @rdname p_lm
#' @export
gen_lm <- function(n, k, R2, R2_0 = 0, k.R2_0 = 0,
				   R2.resid=1-R2, fixed.X=TRUE, ...){
	if(!fixed.X){
		X <- matrix(rnorm(k*n), n, k)
	} else {
		lst <- vector('list', k)
		for(i in 1:k) lst[[i]] <- 0:1
		x <- expand.grid(lst)
		X <- if(nrow(x) < n)
			x[rep(1:nrow(x), each=ceiling(n/nrow(x))), ]
		else
			x[floor(seq(1, nrow(x), length.out=n)),]
		X <- X[1:n, ]
		X <- scale(X)
	}
	colnames(X) <- paste0('X', 1:k)
	R2s <- R2 - R2_0
	betas <- c(sqrt(R2s), sqrt(R2_0), rep(0, k-2))
	y <- colSums(betas * t(X)) + rnorm(n, 0, sqrt(R2.resid))
	data.frame(y, X)
}

# example should include Poisson distribution to match
p_glm <- function(n){

}

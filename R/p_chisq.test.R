#' p-value from chi-squared test simulation
#'
#' Generates multinomial data suitable for analysis with
#' \code{\link{chisq.test}}.
#'
#' @param n sample size per group
#' @param w Cohen's w effect size
#' @param df degrees of freedom
#' @param correct logical; apply continuity correction?
#' @param P0 specific null pattern, specified as a numeric vector or matrix
#' @param P specific power configuration, specified as a numeric vector or matrix
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{matrix} with k rows and k columns
#'   of counts. Default uses \code{\link{gen_chisq.test}}.
#'   User defined version of this function must include the argument \code{...}
#'
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#'
#' @seealso \code{\link{gen_chisq.test}}
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @examples
#'
#' # effect size w + df
#' p_chisq.test(100, w=.2, df=3)
#'
#' # return analysis model
#' p_chisq.test(100, w=.2, df=3, return_analysis=TRUE)
#'
#' # vector of explicit probabilities (goodness of fit test)
#' p_chisq.test(100, P0 = c(.25, .25, .25, .25),
#'                    P = c(.6, .2, .1, .1))
#'
#' # matrix of explicit probabilities (two-dimensional test of independence)
#' p_chisq.test(100, P0 = matrix(c(.25, .25, .25, .25), 2, 2),
#'                    P = matrix(c(.6, .2, .1, .1),2,2))
#'
#' \donttest{
#'     # compare simulated results to pwr package
#'
#'     P0 <- c(1/3, 1/3, 1/3)
#'     P <- c(.5, .25, .25)
#'     w <- pwr::ES.w1(P0, P)
#'     df <- 3-1
#'     pwr::pwr.chisq.test(w=w, df=df, N=100, sig.level=0.05)
#'
#'     # slightly less power when evaluated empirically
#'     p_chisq.test(n=100, w=w, df=df) |> Spower(replications=100000)
#'     p_chisq.test(n=100, P0=P0, P=P) |> Spower(replications=100000)
#'
#'     # slightly differ (latter more conservative due to finite sampling behaviour)
#'     pwr::pwr.chisq.test(w=w, df=df, power=.8, sig.level=0.05)
#'     p_chisq.test(n=interval(50, 200), w=w, df=df) |> Spower(power=.80)
#'     p_chisq.test(n=interval(50, 200), w=w, df=df, correct=FALSE) |>
#'       Spower(power=.80)
#'
#'     # Spower slightly more conservative even with larger N
#'     pwr::pwr.chisq.test(w=.1, df=df, power=.95, sig.level=0.05)
#'     p_chisq.test(n=interval(1000, 2000), w=.1, df=df) |> Spower(power=.95)
#'     p_chisq.test(n=interval(1000, 2000), w=.1, df=df, correct=FALSE) |>
#'            Spower(power=.95)
#'
#' }
#'
#' @export
p_chisq.test <- function(n, w, df, correct = TRUE, P0 = NULL, P = NULL,
						 gen_fun=gen_chisq.test, return_analysis = FALSE, ...) {
	stopifnot(length(n) == 1)
	ret <- if(!missing(w)){
		stopifnot(length(w) == 1)
		stopifnot(length(df) == 1)
		w2 <- w^2
		P0 <- rep(1/(df+1), df+1)
		fn <- function(p1, p0, df, w2){
			ps <- c(p1, rep((1 - p1)/ df, df))
			(sum((ps - p0)^2 / p0) - w2)^2
		}
		# strange that optimize(fn, c(0,1)) gives right w2 but wrong p?
		opt <- optimize(fn, c(P0[1],1), p0=P0, df=df, w2=w2)
		P <- with(opt, c(minimum, rep((1 - minimum)/df, df)))
		# sum((P - P0)^2 / P0) # == w2
		tab <- gen_fun(n=n, P=P, ...)
		chisq.test(tab, correct=correct, p=P0)
	} else {
		tab <- gen_fun(n=n, P=P, ...)
		chisq.test(tab, correct=correct, p=P0)
	}
	if(return_analysis) return(ret)
	p <- ret$p.value
	p
}

#' @rdname p_chisq.test
#' @export
gen_chisq.test <- function(n, P, ...) {
	tab <- as.vector(rmultinom(1, size = n, prob = as.vector(P)))
	if(is.matrix(P))
		tab <- matrix(tab, nrow=nrow(P), ncol=ncol(P))
	tab
}

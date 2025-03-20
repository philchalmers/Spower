#' p-value from proportion test simulation
#'
#' Generates single and multi-sample data
#' for proportion tests and return a p-value. Uses \code{\link{binom.test}}
#' for one-sample applications and \code{\link{prop.test}} otherwise.
#'
#' @param n sample size per group
#' @param h Cohen's h effect size; only supported for one-sample analysis.
#'
#'   Note that it's important to specify the null
#'   value \code{pi} when supplying this effect size as the power
#'   changes depending on these specific values (see example below).
#' @param prop sample probability/proportions of success.
#'   If a vector with two-values or more elements are supplied then
#'   a multi-samples test will be used. Matrices are also supported
#' @param n.ratios allocation ratios reflecting the sample size ratios.
#'   Default of 1 sets the groups to be the same size (n * n.ratio)
#' @param pi probability of success to test against (default is .5). Ignored
#'   for two-sample tests
#' @param two.tailed logical; should a two-tailed or one-tailed test be used?
#' @param exact logical; use fisher's exact test via \code{\link{fisher.test}}?
#'   Use of this flag requires that \code{prop} was specified as a matrix
#' @param correct logical; use Yates' continuity correction?
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{matrix} with two rows and 1 or more
#'   columns. Default uses \code{\link{gen_prop.test}}.
#'   User defined version of this function must include the argument \code{...}
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#'
#' @seealso \code{\link{gen_prop.test}}
#' @return a single p-value
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @examples
#'
#' # one sample, 50 observations, tested against pi = .5 by default
#' p_prop.test(50, prop=.65)
#'
#' # specified using h and pi
#' h <- pwr::ES.h(.65, .4)
#' p_prop.test(50, h=h, pi=.4)
#' p_prop.test(50, h=-h, pi=.65)
#'
#' # two-sample test
#' p_prop.test(50, prop=c(.5, .65))
#'
#' # two-sample test, unequal ns
#' p_prop.test(50, prop=c(.5, .65), n.ratios = c(1,2))
#'
#' # three-sample test, group2 twice as large as others
#' p_prop.test(50, prop=c(.5, .65, .7), n.ratios=c(1,2,1))
#'
#' # Fisher exact test
#' p_prop.test(50, prop=matrix(c(.5, .65, .7, .5), 2, 2))
#'
#' \donttest{
#'     # compare simulated results to pwr package
#'
#'     # one-sample tests
#'     (h <- pwr::ES.h(0.5, 0.4))
#'     pwr::pwr.p.test(h=h, n=60)
#'
#'     # uses binom.test (need to specify null location as this matters!)
#'     Spower(p_prop.test(n=60, h=h, pi=.4))
#'     Spower(p_prop.test(n=60, prop=.5, pi=.4))
#'
#'     # compare with switched null
#'     Spower(p_prop.test(n=60, h=h, pi=.5))
#'     Spower(p_prop.test(n=60, prop=.4, pi=.5))
#'
#'     # two-sample test, one-tailed
#'     (h <- pwr::ES.h(0.67, 0.5))
#'     pwr::pwr.2p.test(h=h, n=80, alternative="greater")
#'     p_prop.test(n=80, prop=c(.67, .5), two.tailed=FALSE,
#'       correct=FALSE) |> Spower()
#'
#'     # same as above, but with continuity correction (default)
#'     p_prop.test(n=80, prop=c(.67, .5), two.tailed=FALSE) |>
#'       Spower()
#'
#'     # three-sample joint test, equal n's
#'     p_prop.test(n=50, prop=c(.6,.4,.7)) |> Spower()
#'
#' }
#'
#' @export
p_prop.test <- function(n, h, prop = NULL, pi = .5,
						n.ratios = rep(1, length(prop)),
						two.tailed = TRUE, correct=TRUE, exact=FALSE,
						gen_fun=gen_prop.test, ...) {
	dat <- gen_prop.test(n=n, h=h, pi=pi, n.ratios=n.ratios, prop=prop, ...)
	if(length(prop) > 1){
		A <- dat[1,]
		B <- dat[2,]
		p <- if(exact){
			stopifnot(is.matrix(prop))
			A <- matrix(A, nrow(prop), ncol(prop))
			fisher.test(A)$p.value
		} else prop.test(A, B, correct=correct)$p.value
	} else {
		p <- binom.test(dat[1,1], n=n, p=pi)$p.value
	}
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' @rdname p_prop.test
#' @export
gen_prop.test <- function(n, h, prop = NULL, pi = .5, n.ratios = rep(1, length(prop)), ...) {
	if(!missing(h)){
		root.h <- function(p1, p2, h)
			(2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))) - h
		n.ratios <- 1
		int <- if(h > 0) c(pi, 1) else c(0, pi)
		root <- uniroot(root.h, interval=int,
						p2=ifelse(length(n.ratios) == 2, prop, pi), h=h)$root
		prop <- root
		# pwr::ES.h(prop, pi) == h
	}
	stopifnot(length(n) == 1)
	n.each <- n * n.ratios
	stopifnot(all.equal(n.each, as.integer(n.each)))
	dat <- if(length(prop) > 1){
		sapply(1:length(prop), \(i){
			vals <- rbinom(n * n.ratios[i], 1, prob=prop[i])
			c(sum(vals), length(vals))
		})
	} else {
		matrix(c(sum(rbinom(n, 1, prob = prop)), n), nrow=2)
	}
	dat
}

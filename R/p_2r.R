#' p-value from comparing two or more correlations simulation
#'
#' Function utilizes \code{\link[cocor]{cocor}} to perform correlation
#' comparison for independent, overlapping, and non-overlapping correlation designs.
#' Type type of correlation design is inferred based on which correlations are
#' specified.
#'
#' For independent group tests, only \code{r.ab} and \code{r.ab2} need to be specified,
#' where the null hypothesis pertains to \eqn{H_0: r_{ab}=r_{ab2}}.
#'
#' For overlapping correlation tests, \code{r.ab}, \code{r.ac}, and \code{r.bc}
#' need to be specified, where the null hypothesis pertains to \eqn{H_0: r_{ab}=r_{ac}}.
#'
#' For non-overlapping correlation tests, all correlations expect for \code{r.ab2} must be
#' specified, where the null hypothesis pertains to \eqn{H_0: r_{ab}=r_{cd}}.
#'
#' @param n sample size
#' @param r.ab correlation between variable A and B (for independent groups,
#'   this is for sample 1)
#' @param r.ab2 (for independent group test only) correlation between
#'   variable A and B in sample 2
#' @param r.ac (for overlap/non-overlap) correlation between A and C. This is the correlation
#'   used in the overlapping hypothesis test, comparing this input to \code{r.ab}
#' @param r.bc (for overlap/non-overlap only) correlation between B and C.
#' @param r.ad (for non-overlap only) correlation between A and D
#' @param r.bd (for non-overlap only) correlation between B and D
#' @param r.cd (for non-overlap only) correlation between C and D. This is the correlation
#'   used in the non-overlapping hypothesis test, comparing this input to \code{r.ab}
#' @param n2_n1 sample size ratio. Only used for independent group test
#' @param two.tailed logical; use two-tailed test?
#' @param test hypothesis testing method to use. Defaults to \code{'fisher1925'}
#'   for the independent groups test and \code{'steiger1980'} for overlap/non-overlap tests
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{matrix} with \code{n} rows.
#'   Default uses \code{\link{gen_2r}}.
#'   User defined version of this function must include the argument \code{...}
#' @param return_analysis logical; return the analysis object for further
#'   extraction and customization?
#' @param ... additional arguments to be passed to \code{gen_fun}. Not used
#'   unless a customized \code{gen_fun} is defined
#' @importFrom cocor cocor
#' @importFrom methods slot
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @return a single p-value
#' @export
#' @examples
#'
#' # independent (same x-y pairing across groups)
#' p_2r(100, r.ab=.5, r.ab2=.6)
#'
#' # return cocor object for further analysis
#' p_2r(100, r.ab=.5, r.ab2=.6, return_analysis = TRUE)
#'
#' \donttest{
#'
#'    # estimate empirical power
#'    p_2r(n=100, r.ab=.5, r.ab2=.6) |> Spower()
#'
#'    # estimate n required to reach 80% power
#'    p_2r(n=interval(100, 5000), r.ab=.5, r.ab2=.6) |> Spower(power=.80)
#'
#' }
#'
#' # overlap (same y, different xs)
#' # H0: r.ab = r.bc
#' p_2r(100, r.ab=.5, r.ac=.3, r.bc=.2)
#'
#' # nonoverlap (different ys, different xs)
#' # H0: r.ab = r.cd
#' p_2r(100, r.ab=.5, r.ac=.3, r.bc=.2, r.ad=.2, r.bd=.4, r.cd=.2)
#'
#'
p_2r <- function(n, r.ab, r.ab2 = NULL, r.ac, r.bc,
				 r.ad, r.bd, r.cd, n2_n1 = 1, two.tailed=TRUE,
				 test=NULL, gen_fun=gen_2r,
				 return_analysis = FALSE, ...){
	type <- 'overlap'
	if(!is.null(r.ab2))
		type <- 'independent'
	if(!missing(r.bd))
		type <- 'nonoverlap'
	if(type == 'independent'){
		if(is.null(test)) test <- 'fisher1925'
		R1 <- matrix(c(1,r.ab, r.ab, 1), 2, 2)
		R2 <- matrix(c(1,r.ab2, r.ab2, 1), 2, 2)
		cnms <- c('y', 'x')
		df1 <- data.frame(gen_fun(n, R=R1, ...))
		df2 <- data.frame(gen_fun(n * n2_n1, R=R2, ...))
		colnames(df1) <- colnames(df2) <- cnms
		dat <- list(sample1=df1, sample2=df2)
	} else if(type == 'overlap'){
		if(is.null(test)) test <- 'steiger1980'
		R <- matrix(c(1,r.ab, r.ac,
					   r.ab, 1, r.bc,
					   r.ac,r.bc, 1), 3, 3)
		dat <- data.frame(gen_fun(n, R=R, ...))
		colnames(dat) <- c('y', 'x1', 'x2')
	} else {
		if(is.null(test)) test <- 'steiger1980'
		R <- matrix(c(1, r.ab, r.ac, r.ad,
					   r.ab, 1, r.bc, r.bd,
					   r.ac, r.bc, 1, r.cd,
					   r.ad, r.bd, r.cd, 1), 4, 4)
		dat <- data.frame(gen_fun(n, R=R, ...))
		colnames(dat) <- c('y1', 'x1', 'y2', 'x2')
	}
	res <- if(type == 'independent'){
		cocor::cocor(~ y + x | y + x, dat, test=test)
	} else if(type == 'overlap'){
		cocor::cocor(~ y + x1 | y + x2, dat, test=test)
	} else {
		cocor::cocor(~ y1 + x1 | y2 + x2, dat, test=test)
	}
	if(return_analysis) return(res)
	pick <- methods::slot(res, test)
	p <- pick$p.value
	p <- ifelse(two.tailed, p, p/2)
	p
}

#' @rdname p_2r
#' @param R a correlation matrix constructed from the inputs
#'   to \code{\link{p_2r}}
#' @export
gen_2r <- function(n, R, ...){
	dat <- SimDesign::rmvnorm(n, mean=numeric(nrow(R)), sigma=R)
	dat
}

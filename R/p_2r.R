#' p-value from comparing two or more correlations simulation
#'
#' Function utilizes \code{\link[cocor]{cocor}} to perform correlation
#' comparison for independent, overlapping, and non-overlapping designs.
#'
#' @param n sample size
#' @param r.ab1 correlation between variable A and B in sample 1
#' @param r.ab2 correlation between variable A and B in sample 2
#' @param r.ac1 same pattern as \code{r.ab1}
#' @param r.ac2 same pattern as \code{r.ab2}
#' @param r.bc1 ...
#' @param r.bc2 ...
#' @param r.ad1 ...
#' @param r.ad2 ...
#' @param r.bd1 ...
#' @param r.bd2 ...
#' @param r.cd1 ...
#' @param r.cd2 ...
#' @param n2_n1 sample size ratio
#' @param two.tailed logical; use two-tailed test?
#' @param type type of correlation design
#' @param test hypothesis method to use. Defaults to 'fisher1925'
#' @param gen_fun function used to generate the required discrete data.
#'   Object returned must be a \code{matrix} with \code{n} rows.
#'   Default uses \code{\link{gen_2r}}.
#'   User defined version of this function must include the argument \code{...}
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
#' p_2r(100, r.ab1=.5, r.ab2=.6)
#'
#' \donttest{
#'
#'    # estimate empirical power
#'    p_2r(n=100, r.ab1=.5, r.ab2=.6) |> Spower()
#'
#'    # estimate n required to reach 80% power
#'    p_2r(n=NA, r.ab1=.5, r.ab2=.6) |>
#'         Spower(power=.80, interval=c(100, 5000))
#'
#' }
#'
#' # overlap (same y, different xs)
#' p_2r(100, r.ab1=.5, r.ab2=.7,
#'           r.ac1=.3, r.ac2=.3,
#'           r.bc1=.2, r.bc2=.2, type = 'overlap')
#'
#' # nonoverlap (different ys, different xs)
#' p_2r(100, r.ab1=.5, r.ab2=.6,
#'           r.ac1=.3, r.ac2=.3,
#'           r.bc1=.2, r.bc2=.2,
#'           r.ad1=.2, r.ad2=.2,
#'           r.bd1=.4, r.bd2=.4,
#'           r.cd1=.2, r.cd2=.2,
#'           type = 'nonoverlap')
#'
#'
p_2r <- function(n, r.ab1, r.ab2, r.ac1, r.ac2, r.bc1, r.bc2,
				 r.ad1, r.ad2, r.bd1, r.bd2, r.cd1, r.cd2,
				 n2_n1 = 1, two.tailed=TRUE,
				 type = c('independent', 'overlap', 'nonoverlap'),
				 test = 'fisher1925', gen_fun=gen_2r, ...){
	type <- match.arg(type)
	if(type == 'independent'){
		R1 <- matrix(c(1,r.ab1, r.ab1, 1), 2, 2)
		R2 <- matrix(c(1,r.ab2, r.ab2, 1), 2, 2)
		cnms <- c('y', 'x')
	} else if(type == 'overlap'){
		R1 <- matrix(c(1,r.ab1, r.ac1,
					   r.ab1, 1, r.bc1,
					   r.ac1,r.bc1, 1), 3, 3)
		R2 <- matrix(c(1,r.ab2, r.ac2,
					   r.ab2, 1, r.bc2,
					   r.ac2, r.bc2, 1), 3, 3)
		cnms <- c('y', 'x1', 'x2')
	} else {
		R1 <- matrix(c(1,r.ab1, r.ac1, r.ad1,
					   r.ab1, 1, r.bc1, r.bd1,
					   r.ac1, r.bc1, 1, r.cd1,
					   r.ad1, r.bd1, r.cd1, 1), 4, 4)
		R2 <- matrix(c(1,r.ab2, r.ac2, r.ad2,
					   r.ab2, 1, r.bc2, r.bd2,
					   r.ac2, r.bc2, 1, r.cd2,
					   r.ad2, r.bd2, r.cd2, 1), 4, 4)
		cnms <- c('y1', 'x1', 'y2', 'x2')
	}
	df1 <- data.frame(gen_fun(n, R=R1, ...))
	df2 <- data.frame(gen_fun(n * n2_n1, R=R2, ...))
	colnames(df1) <- colnames(df2) <- cnms
	dat <- list(sample1=df1, sample2=df2)
	res <- if(type == 'independent'){
		cocor::cocor(~ y + x | y + x, dat, test=test)
	} else if(type == 'overlap'){
		cocor::cocor(~ y + x1 | y + x2, dat, test=test)
	} else {
		cocor::cocor(~ y1 + x1 | y2 + x2, dat, test=test)
	}
	pick <- methods::slot(res, test)
	p <- pick$p.value
	p <- ifelse(two.tailed, p, p*2)
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

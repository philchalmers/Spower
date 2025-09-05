#' Spower evaluations across different combinations
#'
#' The function \code{\link{SpowerBatch}}, on the other hand, can be used to
#' run \code{\link{Spower}} across
#' different simulation combinations, returning a \code{list} of results instead.
#' Can also be used as a pre-computing step before using
#' \code{\link{SpowerCurve}}, and shares the same syntax specification (see
#' \code{\link{SpowerCurve}} for further examples).
#'
#' @rdname Spower
#' @export
#'
#' @examples
#'
#' ##############################################
#' # SpowerBatch() examples
#' ##############################################
#'
#' \dontrun{
#'
#' # estimate power given varying sample sizes
#' p_t.test(d=0.2) |>
#'   SpowerBatch(n=c(30, 90, 270, 550), replications=1000) -> nbatch
#' nbatch
#'
#' # can be stacked to view the output as data.frame
#' as.data.frame(nbatch)
#'
#' # plot with SpowerCurve()
#' SpowerCurve(batch=nbatch)
#'
#' # equivalent, but re-runs the computations
#' p_t.test(d=0.2) |> SpowerCurve(n=c(30, 90, 270, 550), replications=1000)
#'
#' # estimate power given varying sample sizes and effect size
#' p_t.test() |> SpowerBatch(n=c(30, 90, 270, 550),
#'                           d=c(.2, .5, .8), replications=1000) -> ndbatch
#' ndbatch
#'
#' # plot with SpowerCurve()
#' SpowerCurve(batch=ndbatch)
#'
#' }
#'
#'
SpowerBatch <- function(..., interval = NULL, power = NA,
						sig.level=.05, beta_alpha = NULL, sig.direction = 'below',
						replications=10000, integer,
						parallel = FALSE, cl = NULL,
						ncores = parallelly::availableCores(omit = 1L),
						predCI = 0.95, predCI.tol = .01, verbose = TRUE,
						check.interval=FALSE, maxiter=150, wait.time = NULL,
						select = NULL, control = list()){

	dots <- match.call(expand.dots = FALSE)$...
	if(all(is.na(sig.level))){
		interval <- c(0,1)
		integer <- FALSE
	}
	expr <- dots[[1]]
	expr <- match.call(eval(expr[[1]], envir = parent.frame(1)), expr)
	pick <- if(length(dots) > 1) names(dots[-1]) else NULL
	if(all(is.na(power))){
		conditions <- do.call(SimDesign::createDesign, c(dots[-1], sig.level=sig.level, power=power))
	} else {
		if(is.null(interval))
			stop('search interval must be included', call.=FALSE)
		lst_expr <- as.list(expr)[-1]
		if(length(lst_expr))
			lst_expr <- lst_expr[sapply(lst_expr, \(x) is.atomic(x) || is.list(x))]
		conditions <- do.call(SimDesign::createDesign, c(lst_expr,
														 dots[-1],
														 sig.level=list(sig.level),
														 power=list(power)))
	}
	if(!all(rowSums(is.na(conditions)) == 1))
		stop('Exactly *one* argument must be set to \'NA\' in SpowerBatch(..., power, sig.level)',
			 call.=FALSE)
	power <- conditions$power
	sig.level <- conditions$sig.level
	if(!is.na(power[1])){
		if(missing(integer)){
			integer <- !(has.decimals(interval) || diff(interval) < 5)
			if(!integer && verbose)
				message('\nUsing continuous search interval (integer = FALSE).')
		}
	} else integer <- FALSE
	control$nparent <- 2
	out <- vector('list', nrow(conditions))
	names(out) <- paste0('CONDITION_', 1:length(out))
	for(i in 1:length(out)){
		row <- conditions[i, ]
		tmpexpr <- expr
		if(length(pick))
			tmpexpr[pick] <- row[,pick]
		out[[i]] <- do.call(Spower, c(tmpexpr,
									  list(power=power[i], sig.level=sig.level[i], beta_alpha=beta_alpha,
									  	 interval=interval, integer=integer, replications=replications,
									  	 parallel=parallel, cl=cl, predCI=predCI, predCI.tol=predCI.tol,
									  	 verbose=verbose, check.interval=check.interval, sig.direction=sig.direction,
									  	 maxiter=maxiter, wait.time=wait.time, select=select, control=control)))
	}
	class(out) <- 'SpowerBatch'
	attr(out[[1]], 'Spower_extra')$full_conditions <- conditions
	.SpowerEnv$lastSim <- out
	out
}

#' @rdname Spower
#' @export
print.SpowerBatch <- function(x, ...){
	class(x) <- 'list'
	print(x, ...)
}

#' @rdname Spower
#' @export
as.data.frame.SpowerBatch <- function(x, ...){
	y <- lapply(x, as.data.frame)
	out <- do.call(rbind, y)
	rownames(out) <- NULL
	out
}




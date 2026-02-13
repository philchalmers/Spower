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
#'
#' @param fully.crossed logical; should the supplied conditions to
#'   \code{SpowerBatch} be fully crossed? Passed to the same
#'   argument documented in \code{\link[SimDesign]{createDesign}}
#'
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
#' # For non-crossed experimental combinations, pass fully.crossed = FALSE. Note
#' # that this requires the lengths of the inputs to match
#' p_t.test() |> SpowerBatch(n=c(30, 90, 270),
#'            d=c(.2, .5, .8), replications=1000, fully.crossed=FALSE) -> batch3
#'
#'
#' ##############################
#'
#' # Batches also useful for drawing graphics outside of current framework
#' # in SpowerCurve(). Below an image is drawn pertaining to the distribution
#' # of the effects (H0 vs Ha hypotheses), giving the classic sampling
#' # distribution comparisons of the effect sizes, however presents the
#' # information using kernel density plots as this may be useful when the
#' # sampling distributions are non-normal
#'
#' # Define wrapper function that returns p-value and estimated mean difference
#' Ice_T <- function(...){
#' 	  out <- p_t.test(..., return_analysis=TRUE)
#' 	  ret <- c(p=out$p.value, mu_d=unname(with(out, estimate[1] - estimate[2])))
#'    ret
#' }
#'
#' # rapper returns p-value and effect size of interest
#' Ice_T(n=90, d=.5)
#'
#' # run batch mode to get 4 mean difference combinations, selecting out only
#' # the 'p' for the power-analysis portions
#' batch <- Ice_T(n=90) |>
#'    SpowerBatch(d=c(0, .2, .5, .8), select="p")
#' batch
#' as.data.frame(batch)
#'
#' # create big table of results across the batches
#' results <- SimResults(batch, rbind=TRUE)
#' results$d <- factor(results$d)
#' results
#'
#' # draw H0 vs Ha relationships for each effect size
#' library(ggplot2)
#' library(patchwork)
#' gg1 <- ggplot(subset(results, d %in% c(0, .2)),
#' 			  aes(mu_d, colour=d)) +
#' 	geom_density() + ggtitle('Small effect (d = 0.2)') +
#' 	theme(legend.position='none') +
#' 	xlab(expression(mu[d])) + xlim(c(-0.75, 1.5))
#' gg2 <- ggplot(subset(results, d %in% c(0, .5)),
#' 			  aes(mu_d, colour=d)) +
#' 	  geom_density() + ggtitle('Medium effect  (d = 0.5)') +
#' 	  theme(legend.position='none') + xlab(expression(mu[d])) +
#' 	  xlim(c(-0.75, 1.5))
#' gg3 <- ggplot(subset(results, d %in% c(0, .8)),
#' 			  aes(mu_d, colour=d)) +
#' 	  geom_density() + ggtitle('Large effect  (d = 0.8)') +
#' 	  theme(legend.position='none') + xlab(expression(mu[d])) +
#' 	  xlim(c(-0.75, 1.5))
#'
#' gg1 / gg2 / gg3
#'
#' }
#'
#'
SpowerBatch <- function(..., interval = NULL, power = NA,
						sig.level=.05, beta_alpha = NULL, sig.direction = 'below',
						replications=10000, integer, fully.crossed = TRUE,
						parallel = FALSE, cl = NULL,
						ncores = parallelly::availableCores(omit = 1L),
						predCI = 0.95, predCI.tol = .01, verbose = interactive(),
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
		conditions <- do.call(SimDesign::createDesign, c(dots[-1], sig.level=sig.level, power=power,
														 fully.crossed=fully.crossed))
	} else {
		if(is.null(interval))
			stop('search interval must be included', call.=FALSE)
		lst_expr <- as.list(expr)[-1]
		if(length(lst_expr))
			lst_expr <- lst_expr[sapply(lst_expr, \(x) is.atomic(x) || is.list(x))]
		conditions <- do.call(SimDesign::createDesign, c(lst_expr,
														 dots[-1],
														 sig.level=list(sig.level),
														 power=list(power),
														 fully.crossed=fully.crossed))
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




#' Draw power curve from simulation functions
#'
#' Draws power curves that either a) estimate the power given a
#' set of varying conditions or b) solves a set of root conditions
#' given fixed values of power. Confidence/prediction intervals are
#' included in the output to reflect the estimate uncertainties.
#'
#' @param sim function that both creates the data and returns a single
#'   p-value for the analysis of interest
#'
#' @param ... a set of conditions to use in the simulation that must match the
#'   arguments in the function \code{sim}. See \code{\link{Spower}}
#'
#' @param power power level to use. If set to \code{NA} then the empirical power
#'   will be estimated given the fixed \code{...} input; otherwise,
#'   can be specified as a vector to solve the missing elements in
#'   \code{...}
#'
#' @param varying either a vector of values to substitute into the missing \code{...}
#'   terms or a structure created from \code{\link[SimDesign]{createDesign}}.
#'   The benefit of the ladder is that multiple factors can vary, and NA
#'   placeholders will not be required. Requires \code{power} to be set to \code{NA}.
#'
#'   Note that only the first three columns in this object will be plotted using
#'   the x-y, colour, and facet wrap aesthetics, respectively. However,
#'   if necessary the data can be extracted for further visualizations via
#'   \code{\link[ggplot2]{ggplot_build}} to provide more customized control
#'
#' @param maxiter see \code{\link{Spower}}
#'
#' @param sig.level see \code{\link{Spower}}
#'
#' @param interval search interval to use when \code{\link{SimSolve}} is required.
#'   Can be a vector of length two to apply the same interval across
#'   the \code{varying} information or a \code{matrix} with two columns
#'   to apply intervals on a per-row basis
#'
#' @param plotCI logical; include confidence/prediction intervals in plots?
#'
#' @param wait.time see \code{\link{Spower}}
#' @param replications see \code{\link{Spower}}
#' @param integer see \code{\link{Spower}}
#' @param parallel see \code{\link{Spower}}
#' @param cl see \code{\link{Spower}}
#' @param ncores see \code{\link{Spower}}
#' @param control see \code{\link{Spower}}
#' @param predCI see \code{\link{Spower}}
#' @param predCI.tol see \code{\link{Spower}}
#' @param check.interval see \code{\link{Spower}}, though is set to \code{FALSE}
#'   by default instead
#' @param verbose see \code{\link{Spower}}
#' @param prior see \code{\link{Spower}}
#' @return a ggplot2 object
#' @import ggplot2
#' @export
#'
#' @seealso \code{\link{Spower}}
#'
#' @examples
#' \dontrun{
#'
#' # estimate power given varying sample sizes
#' gg <- powerCurve(p_t.test, varying=c(30, 90, 270, 550), n=NA, d=0.2,
#'  		   replications=1000)
#'
#' # alternatively, specifying varying as a createDesign() object (does not
#' # require NA placeholders)
#' varying <- createDesign(n=c(30, 90, 270, 550))
#' gg <- powerCurve(p_t.test, varying=varying, d=0.2, replications=1000)
#'
#' # also equivalent, though no CIs plotted
#' varying <- createDesign(n=c(30, 90, 270, 550),
#'                         d=0.2)
#' gg <- powerCurve(p_t.test, varying=varying, replications=1000, plotCI=FALSE)
#'
#' #####
#'
#' # Because output is a ggplot2 object can be modified
#' gg + geom_text(aes(label=power), size=5, colour='red', nudge_y=.05) +
#'   ylab(expression(1-beta)) + theme_grey()
#'
#' # using default precision (10000 replications). Parallel computations
#' # generally recommended in this case to save time
#' powerCurve(p_t.test, varying=c(30, 90, 270, 550), n=NA, d=0.2)
#'
#' # alternatively, get information from last printed graphic instead of saving
#' gg <- last_plot()
#' gg + coord_flip() # flip coordinates
#'
#' # estimate sample sizes given varying power
#' powerCurve(p_t.test, n=NA, d=0.2, interval=c(10, 1000),
#' 		   power=c(.1, .25, .5, .75, .9), maxiter=30)
#'
#' # estimate power varying d
#' powerCurve(p_t.test, varying=seq(.1, 1, by=.2), n=50, d=NA,
#' 		   replications=1000)
#'
#' # estimate d varying power
#' powerCurve(p_t.test, n=50, d=NA,
#' 		   maxiter=30, interval=c(.01, 1),
#' 		   power=c(.1, .25, .5, .75, .9))
#'
#' #####
#'
#' # vary two inputs instead of one (second column uses colour aesthetic)
#' varying <- createDesign(n=c(30, 90, 270, 550),
#'                         d=c(.2, .5, .8))
#' powerCurve(p_t.test, varying=varying, replications=2000)
#'
#' # extract data for alternative presentations
#' build <- ggplot_build(last_plot())
#' build
#'
#' df <- build$plot$data
#' head(df)
#' ggplot(df, aes(n, power, linetype=d)) + geom_line()
#'
#' # vary three inputs (third column uses facet_wrap)
#' varying <- createDesign(n=c(30, 90, 270),
#'                         var.equal=c(FALSE, TRUE),
#'                         d=c(.2, .5))
#' powerCurve(p_t.test, varying=varying, replications=1000, plotCI=FALSE)
#'
#' }
#'
powerCurve <- function(sim, varying, ..., interval = NULL, power = NA,
					   sig.level=.05, replications=10000, integer, prior = NULL,
					   plotCI=TRUE, parallel = FALSE, cl = NULL,
					   ncores = parallelly::availableCores(omit = 1L),
					   predCI = 0.95, predCI.tol = .01, verbose = TRUE,
					   check.interval=FALSE, maxiter=150, wait.time = NULL,
					   control = list()){
	dots <- dotse <- list(...)
	opower <- power
	if(!missing(varying) && is.numeric(varying)){
		pick <- sapply(dots, \(x) all(is.na(x)))
		column <- names(dots)[pick]
		varying <- data.frame(varying)
		colnames(varying) <- column
	}
	if(!missing(varying))
		column <- colnames(varying)
	if(is.na(sig.level))
		stop('NA for sig.level not supported')
	if(length(power) > 1 && !missing(varying))
		stop('Either power is fixed to NA or varying is specified, not both')
	if(all(is.na(power)) && missing(varying))
		stop('Must specify varying')
	if(length(opower) == 1)
		opower <- rep(opower, nrow(varying))
	out <- vector('list', length(opower))
	for(i in 1:length(out)){
		dotse <- dots
		power <- opower[i]
		if(is.na(power))
			dotse[column] <- varying[i,]
		if(length(power) == 1 && !is.na(power) && !is.na(sig.level)){
			if(missing(integer)){
				integer <- !(has.decimals(interval) || diff(interval) < 5)
				if(!integer && verbose)
					message('\nUsing continuous search interval (set manually by passing integer = FALSE).')
			}
		} else integer <- FALSE
		out[[i]] <- do.call(Spower, c(sim=sim, dotse,
									  list(power=power, sig.level=sig.level, beta_alpha=NULL,
									  interval=interval, integer=integer, replications=replications,
									  parallel=parallel, cl=cl, predCI=predCI, predCI.tol=predCI.tol,
									  verbose=verbose, check.interval=check.interval,
									  maxiter=maxiter, wait.time=wait.time, control=control)))
	}
	CI.low <- CI.high <- NULL # for check?
	main <- "Power Curve"
	if(is.na(power)){
		CI <- unname(t(sapply(out, \(x) summary(x)$power.CI)))
		df <- data.frame(do.call(rbind, out), CI.low=CI[,1], CI.high=CI[,2])
		if(ncol(varying) > 1 && nrow(unique(varying[,2])) > 1){
			df[[column[2]]] <- factor(df[[column[2]]])
			gg <- ggplot(df, aes(.data[[column[1]]], power,
								 color=.data[[column[2]]]))
			if(plotCI){
				gg <- gg + geom_ribbon(aes(ymin=CI.low, ymax=CI.high),
									   alpha=.2, linetype='dashed')
				main <- "Power Curve (with 95% CIs)"
			}
			gg <- gg + geom_line() + geom_point() +	ggtitle(main) +	theme_bw()
			if(ncol(varying) > 2 && nrow(unique(varying[,3])) > 1)
				gg <- gg + facet_wrap( ~ .data[[column[3]]])
		} else {
			gg <- ggplot(df, aes(.data[[column[1]]], power))
			if(plotCI){
				gg <- gg + geom_ribbon(aes(ymin=CI.low, ymax=CI.high), alpha=.2) +
					geom_line(aes(y=CI.low), linetype='dashed') +
					geom_line(aes(y=CI.high), linetype='dashed')
				main <- "Power Curve (with 95% CIs)"
			}
			gg <- gg + geom_line() + geom_point() +	ggtitle(main) +	theme_bw()
		}
	} else {
		CI <- unname(t(sapply(out, \(x) summary(x)$predCIs_root)))
		df <- cbind(do.call(rbind, out), CI.low=CI[,1], CI.high=CI[,2])
		pick <- sapply(dots, \(x) all(is.na(x)))
		column <- names(dots)[pick]
		gg <- ggplot(df, aes(power, .data[[column]]))
		if(plotCI){
			gg <- gg + geom_ribbon(aes(ymin=CI.low, ymax=CI.high), alpha=.2) +
				geom_line(aes(y=CI.low), linetype='dashed') +
				geom_line(aes(y=CI.high), linetype='dashed')
			main <- "Power Curve (with 95% PIs)"
		}
		gg <- gg + geom_line() + geom_point() +	ggtitle(main) +	theme_bw()
	}
	print(gg)
	invisible(gg)
}

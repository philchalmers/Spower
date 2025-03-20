#' Draw power curve from simulation functions
#'
#' Draws power curves that either a) estimate the power given a
#' set of varying conditions or b) solves a set of root conditions
#' given fixed values of power. Confidence/prediction intervals are
#' included in the output to reflect the estimate uncertainties, though note
#' that fewer replications/iterations are used compared to
#' \code{\link{Spower}} as the goal is visualization of competing
#' variable inputs rather than precision of a given input.
#'
#'
#' @param ... first expression input must be identical to \code{...} in
#' \code{\link{Spower}}, while the remaining named inputs must match the arguments
#' to this expression to indicate which variables should be modified in the
#' resulting power curves. Providing \code{NA} values is also supported to
#' solve the missing component
#'
#' Note that only the first three named arguments will be plotted using
#' the x-y, colour, and facet wrap aesthetics, respectively. However,
#' if necessary the data can be extracted for further visualizations via
#' \code{\link[ggplot2]{ggplot_build}} to provide more customized control
#'
#' @param power power level to use. If set to \code{NA} then the empirical power
#'   will be estimated given the fixed \code{...} input; otherwise,
#'   can be specified as a vector to solve the missing elements in
#'   \code{...}
#'
#' @param maxiter see \code{\link{Spower}}, though set to 50 instead of 150
#'
#' @param sig.level see \code{\link{Spower}}
#'
#' @param interval search interval to use when \code{\link[SimDesign]{SimSolve}} is required.
#'   Can be a vector of length two to apply the same interval across
#'   the \code{varying} information or a \code{matrix} with two columns
#'   to apply intervals on a per-row basis
#'
#' @param plotCI logical; include confidence/prediction intervals in plots?
#'
#' @param wait.time see \code{\link{Spower}}
#' @param replications see \code{\link{Spower}}, though set to 2500 instead of
#'   10000
#' @param integer see \code{\link{Spower}}
#' @param parallel see \code{\link{Spower}}
#' @param cl see \code{\link{Spower}}
#' @param ncores see \code{\link{Spower}}
#' @param control see \code{\link{Spower}}
#' @param predCI see \code{\link{Spower}}
#' @param predCI.tol see \code{\link{Spower}}
#' @param check.interval see \code{\link{Spower}}, though is set to \code{FALSE}
#'   by default instead
#' @param plotly logical; draw the graphic into the interactive \code{plotly}
#'   interface? If \code{FALSE} the ggplot2 object will be returned instead
#' @param verbose see \code{\link{Spower}}
#' @return a ggplot2 object
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @export
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#'
#' @seealso \code{\link{Spower}}
#'
#' @examples
#' \donttest{
#'
#' # estimate power given varying sample sizes
#' gg <- p_t.test(d=0.2) |> SpowerCurve(n=c(30, 90, 270, 550))
#'
#' # Output is a ggplot2 (rendered with plotly by default);hence, can be modified
#' gg + geom_text(aes(label=power), size=5, colour='red', nudge_y=.05) +
#'   ylab(expression(1-beta)) + theme_grey()
#'
#' # Increase precision by using 10000 replications. Parallel computations
#' #   generally recommended in this case to save time
#' p_t.test(d=0.2) |> SpowerCurve(n=c(30, 90, 270, 550), replications=10000)
#'
#' # estimate sample sizes given varying power
#' p_t.test(n=NA, d=0.2) |>
#'   SpowerCurve(power=c(.2, .4, .6, .8), interval=c(10, 1000))
#'
#' # get information from last printed graphic instead of saving
#' gg <- last_plot()
#' gg + coord_flip() # flip coordinates to put power on y-axis
#'
#' # estimate power varying d
#' p_t.test(n=50) |> SpowerCurve(d=seq(.1, 1, by=.2))
#'
#' # estimate d varying power
#' p_t.test(n=50, d=NA) |>
#'   SpowerCurve(power=c(.2, .4, .6, .8), interval=c(.01, 1))
#'
#'
#' #####
#'
#' # vary two inputs instead of one (second input uses colour aesthetic)
#' p_t.test() |> SpowerCurve(n=c(30, 90, 270, 550),
#'                          d=c(.2, .5, .8))
#'
#' # extract data for alternative presentations
#' build <- ggplot_build(last_plot())
#' build
#'
#' df <- build$plot$data
#' head(df)
#' ggplot(df, aes(n, power, linetype=d)) + geom_line()
#'
#' # vary three arguments (third uses facet_wrap ... any more than that and
#' #   you're on your own!)
#' SpowerCurve(p_t.test, varying=varying)
#'
#' p_t.test() |> SpowerCurve(n=c(30, 90, 270, 550),
#'                          d=c(.2, .5, .8),
#'                          var.equal=c(FALSE, TRUE))
#'
#' }
#'
SpowerCurve <- function(..., interval = NULL, power = NA,
					   sig.level=.05, replications=2500, integer,
					   plotCI=TRUE, plotly=TRUE, parallel = FALSE, cl = NULL,
					   ncores = parallelly::availableCores(omit = 1L),
					   predCI = 0.95, predCI.tol = .01, verbose = TRUE,
					   check.interval=FALSE, maxiter=50, wait.time = NULL,
					   control = list()){
	dots <- match.call(expand.dots = FALSE)$...
	if(is.na(sig.level))
		stop('solving for sig.level not yet supported', call.=FALSE)
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
		stop('Exactly *one* argument must be set to \'NA\' in SpowerCurve(..., power, sig.level)',
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
	for(i in 1:length(out)){
		row <- conditions[i, ]
		tmpexpr <- expr
		if(length(pick))
			tmpexpr[pick] <- row[,pick]
		out[[i]] <- do.call(Spower, c(tmpexpr,
									  list(power=power[i], sig.level=sig.level[i], beta_alpha=NULL,
									  interval=interval, integer=integer, replications=replications,
									  parallel=parallel, cl=cl, predCI=predCI, predCI.tol=predCI.tol,
									  verbose=verbose, check.interval=check.interval,
									  maxiter=maxiter, wait.time=wait.time, control=control)))
	}
	CI.low <- CI.high <- NULL # for check?
	main <- "Power Curve"
	cnms <- colnames(conditions)
	if(is.na(conditions$power[1])){
		CI <- unname(t(sapply(out, \(x) summary(x)$power.CI)))
		df <- data.frame(do.call(rbind, out), CI.low=CI[,1], CI.high=CI[,2])
		if(length(cnms) > 3){  # more than one varying
			df[[cnms[2]]] <- factor(df[[cnms[2]]])
			gg <- ggplot(df, aes(.data[[cnms[1]]], power,
								 color=.data[[cnms[2]]],
								 fill=.data[[cnms[2]]]))
			if(plotCI){
				gg <- gg + geom_ribbon(aes(ymin=CI.low, ymax=CI.high),
									   alpha=.2, linetype='dashed')
			}
			gg <- gg + geom_line() + geom_point() +
				ggtitle(main) +	theme_bw()
			if(length(cnms) > 4)  # more than two varying
				gg <- gg + facet_wrap( ~ .data[[cnms[3]]])
		} else {   # one varying
			gg <- ggplot(df, aes(.data[[cnms[1]]], power))
			if(plotCI){
				gg <- gg + geom_ribbon(aes(ymin=CI.low, ymax=CI.high), alpha=.2) +
					geom_line(aes(y=CI.low), linetype='dashed') +
					geom_line(aes(y=CI.high), linetype='dashed')
			}
			gg <- gg + geom_line() + geom_point() + ggtitle(main) +
				theme_bw() + theme(plot.title = element_text(hjust = 0.5))
		}
	} else if(var(conditions$power) == 0){
		CI <- unname(t(sapply(out, \(x) summary(x)$predCIs_root)))
		df <- cbind(do.call(rbind, out), CI.low=CI[,1], CI.high=CI[,2])
		pick <- cnms[is.na(conditions[1,])]
		pick2 <- cnms[min(which(pick != cnms))]
		main <- sprintf("Constant Power (%.3f)", power[1])
		if(var(conditions$sig.level) > 0){


		} else {
			gg <- ggplot(df, aes(.data[[pick2]], .data[[pick]]))
			if(plotCI){
				gg <- gg + geom_ribbon(aes(ymin=CI.low, ymax=CI.high), alpha=.2) +
					geom_line(aes(y=CI.low), linetype='dashed') +
					geom_line(aes(y=CI.high), linetype='dashed')
			}
			gg <- gg + geom_line() + geom_point() + ggtitle(main) +
				theme_bw() + theme(plot.title = element_text(hjust = 0.5))
		}
	} else {
		CI <- unname(t(sapply(out, \(x) summary(x)$predCIs_root)))
		df <- cbind(do.call(rbind, out), CI.low=CI[,1], CI.high=CI[,2])
		pick <- cnms[is.na(conditions[1,])]
		gg <- ggplot(df, aes(power, .data[[pick]]))
		if(plotCI){
			gg <- gg + geom_ribbon(aes(ymin=CI.low, ymax=CI.high), alpha=.2) +
				geom_line(aes(y=CI.low), linetype='dashed') +
				geom_line(aes(y=CI.high), linetype='dashed')
		}
		gg <- gg + geom_line() + geom_point() + ggtitle(main) +
			theme_bw() + theme(plot.title = element_text(hjust = 0.5))
	}
	if(plotly){
		print(plotly::ggplotly(gg))
		return(invisible(gg))
	}
	gg
}

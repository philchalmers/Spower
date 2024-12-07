#' @import ggplot2
#' @export
powerCurve <- function(sim, varying, ..., interval = NULL, power = NA,
					   sig.level=.05, replications=10000, integer, prior = NULL,
					   parallel = FALSE, cl = NULL,
					   ncores = parallelly::availableCores(omit = 1L),
					   predCI = 0.95, predCI.tol = .01, verbose = TRUE,
					   check.interval=TRUE, maxiter=150, wait.time = NULL,
					   control = list()){
	dots <- dotse <- list(...)
	opower <- power
	if(is.na(sig.level))
		stop('NA for sig.level not supported')
	if(length(power) > 1 && !missing(varying))
		stop('Either power is fixed to NA or varying is specified, not both')
	if(all(is.na(power)) && missing(varying))
		stop('Must specify varying')
	if(length(opower) == 1)
		opower <- rep(opower, length(varying))
	out <- vector('list', length(opower))
	for(i in 1:length(out)){
		dotse <- dots
		power <- opower[i]
		if(is.na(power)){
			pick <- sapply(dots, \(x) all(is.na(x)))
			dotse[[which(pick)]] <- varying[i]
		}
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
	if(is.na(power)){
		CI <- unname(t(sapply(out, \(x) summary(x)$power.CI)))
		df <- data.frame(do.call(rbind, out), CI.low=CI[,1], CI.high=CI[,2])
		pick <- sapply(dots, \(x) all(is.na(x)))
		column <- names(dots)[pick]
		gg <- ggplot(df, aes(.data[[column]], power)) +
			geom_ribbon(aes(ymin=CI.low, ymax=CI.high), alpha=.3) +
			geom_line() + geom_point() +
			ggtitle("Power Curve (with 95% CIs)")
	} else {
		CI <- unname(t(sapply(out, \(x) summary(x)$predCIs_root)))
		df <- cbind(do.call(rbind, out), CI.low=CI[,1], CI.high=CI[,2])
		pick <- sapply(dots, \(x) all(is.na(x)))
		column <- names(dots)[pick]
		gg <- ggplot(df, aes(power, .data[[column]])) +
			geom_ribbon(aes(ymin=CI.low, ymax=CI.high), alpha=.3) +
			geom_line() + geom_point() +
			ggtitle("Power Curve (with 95% PIs)")
	}
	gg
}

if(FALSE){

	# estimate power given varying sample sizes
	powerCurve(p_t.test, varying=c(30, 60, 90), n=NA, d=0.2,
			   replications=1000)

	# estimate sample sizes given varying power
	powerCurve(p_t.test, n=NA, d=0.2, interval=c(10, 1000),
			   power=c(.1, .25, .5, .75, .9), maxiter=30)

	# estimate power varying d
	powerCurve(p_t.test, varying=seq(.1, 1, by=.1), n=50, d=NA,
			   replications=1000)

	# estimate d varying power
	powerCurve(p_t.test, n=50, d=NA,
			   maxiter=30, interval=c(.01, 1),
			   power=c(.1, .25, .5, .75, .9))

}

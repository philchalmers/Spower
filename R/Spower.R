#' Simulation-based Power Analyses
#'
#' General purpose function that provides power-focused estimates for
#' a priori, prospective/post-hoc, compromise, sensitivity, and criterion power analysis.
#' Function provides a general wrapper to the
#' \code{SimDesign} package's \code{\link[SimDesign]{runSimulation}} and
#' \code{\link[SimDesign]{SimSolve}} functions. As such, parallel processing is
#' automatically supported, along with progress bars,
#' confidence/prediction intervals for the results estimates, safety checks,
#' and more.
#'
#' Five types of power analysis flavors can be performed with \code{Spower},
#' which are triggered based on which supplied input is set to missing (\code{NA}):
#'
#' \describe{
#'    \item{A Priori}{Solve for a missing sample size component
#'      (e.g., \code{n}) to achieve a specific target power rate}
#'    \item{Prospective and Post-hoc}{Estimate the power rate given a set of fixed conditions.
#'      If estimates of effect sizes and other empirical characteristics (e.g., observed sample size)
#'      are supplied this results in observed/retrospective power (not recommended), while if only
#'      sample size is included as the observed quantity, but the effect sizes are treated as unknown, then this
#'      results in post-hoc power (Cohen, 1988)}
#'    \item{Sensitivity}{Solve a missing effect size value as a function of
#'      the other supplied constant components}
#'    \item{Criterion}{Solve the error rate (argument \code{sig.level}) as a
#'      function of the other supplied constant components}
#'    \item{Compromise}{Solve a Type I/Type II error trade-off ratio as a
#'      function of the other supplied constant components and the
#'      target ratio \eqn{q = \beta/\alpha} (argument \code{beta_alpha})}
#' }
#'
#' Prospective and compromise analyses utilize the
#' \code{\link[SimDesign]{runSimulation}} function, while the remaining three
#' approaches utilize the stochastic root solving methods in the function
#' \code{\link[SimDesign]{SimSolve}}.
#' See the example below for a demonstration with an independent samples t-test
#' analysis.
#'
#' @param ... expression to use in the simulation that returns a \code{numeric}
#'   vector containing only p-value information, where the first p-value
#'   in this vector is treated as the focus for all analyses other than prospective/post-hoc power,
#'   or a similarly structure \code{logical} vector when utilizing confidence intervals (CIs).
#'
#'   Internally the first expression is passed to either \code{\link[SimDesign]{SimSolve}} or
#'  \code{\link[SimDesign]{runSimulation}} depending on which element (including
#'  the \code{power} and \code{sig.level} arguments) is set to \code{NA}. For instance,
#'  \code{Spower(p_t.test(n=50, d=.5))} will perform a prospective/post-hoc power evaluation since
#'  \code{power = NA} by default, while \code{Spower(p_t.test(n=NA, d=.5), power = .80)}
#'  will perform an a priori power analysis to solve the missing \code{n} argument.
#'
#'  For expected power computations the arguments to this expression can themselves
#'  be specified as a function to reflect the prior uncertainty. For instance, if
#'  \code{d_prior <- function() rnorm(1, mean=.5, sd=1/8)} then
#'  \code{Spower(p_t.test(n=50, d=d_prior())} will compute the expected power
#'  over the prior sampling distribution for \code{d}
#'
#' @param power power level to use. If set to \code{NA} then the empirical power
#'   will be estimated given the fixed \code{...} inputs
#'   (e.g., for prospective/post-hoc power analysis)
#'
#' @param maxiter maximum number of stochastic root-solving iterations
#'
#' @param select a character vector indicating which elements to
#'   extract from the provided stimulation experiment function. By default, all elements
#'   from the provided function will be used, however if the provided function contains
#'   information not relevant to the power computations (e.g., parameter estimates,
#'   standard errors, etc) then these should be ignored. To extract the complete
#'   results post-analysis use \code{\link[SimDesign]{SimResults}} to allow manual
#'   summarizing of the stored results (applicable only with prospective/post-hoc power)
#'
#' @param sig.level alpha level to use. If set to \code{NA} then the empirical
#'   alpha will be estimated given the fixed \code{conditions} input
#'   (e.g., for criterion power analysis). Only used when the value returned
#'   from the experiment is a \code{numeric} (p-value).
#'
#'   If the return of the supplied experiment is a
#'   \code{logical}, which generally indicates that a confidence interval (CI)
#'   approach were used,
#'   then an argument such as \code{conf.level} should be included
#'   in the simulation experiment to indicate the explicit CI
#'   criteria, and so that this can be modified directly as well
#'
#' @param interval search interval to use when \code{\link[SimDesign]{SimSolve}} is required.
#'   Note that for compromise analyses, where the \code{sig.level} is set to
#'   \code{NA}, if not set explicitly then the interval will default to \code{c(0,1)}
#'
#' @param wait.time (optional) argument to indicate the time to wait
#'  (specified in minutes if supplied as a numeric vector).
#'  See \code{\link[SimDesign]{SimSolve}} for details and
#'  See \code{\link[SimDesign]{timeFormater}} for further specifications
#'
#' @param replications number of replications to use when
#'   \code{\link[SimDesign]{runSimulation}} is required
#'
#' @param lastSpower a previously returned \code{Spower} object to be updated.
#'   Use this if you want to continue where an estimate left off but wish to increase the
#'   precision (e.g., by adding more replications, or by letting the stochastic root solver
#'   continue searching).
#'
#'   Note that if the object was not stored use \code{\link{getLastSpower}}
#'   to obtain the last estimated power object
#'
#' @param integer a logical value indicating whether the search iterations
#'   use integers or doubles.
#'
#'   If missing, automatically set to \code{FALSE} if \code{interval} contains
#'   non-integer numbers or the range is less than 5, as well as
#'   when \code{sig.level = NA}, though in general this should be set explicitly
#'
#' @param beta_alpha (optional) ratio to use in compromise analyses corresponding to
#'   the Type II errors (beta) over the Type I error (alpha). Ratios greater
#'   than 1 indicate that Type I errors are worse than Type II, while ratios
#'   less than one the opposite. A ratio equal to 1 gives an equal trade-off
#'   between Type I and Type II errors
#'
#' @param parallel for parallel computing for slower simulation experiments
#'   (see \code{\link[SimDesign]{runSimulation}} for details).
#'
#   Note that
#   the defined cluster object will be made globally available via
#   \code{\link{SpowerCluster}} so that the cluster definition can be reused
#   automatically
#'
#' @param cl see \code{\link[SimDesign]{runSimulation}}
#'
#' @param ncores see \code{\link[SimDesign]{runSimulation}}
#'
#' @param packages see \code{\link[SimDesign]{runSimulation}}
#'
#' @param control a list of control parameters to pass to
#'   \code{\link[SimDesign]{runSimulation}} or \code{\link[SimDesign]{SimSolve}}
#'
#' @param predCI predicting confidence interval level
#'   (see \code{\link[SimDesign]{SimSolve}})
#'
#' @param predCI.tol predicting confidence interval consistency tolerance
#'    for stochastic root solver convergence (see \code{\link[SimDesign]{SimSolve}}).
#'    Default converges when the power rate CI is consistently
#'    within \code{.01/2} of the target power
#'
#' @param check.interval logical; check the interval range validity
#'   (see \code{\link[SimDesign]{SimSolve}}). Disabled by default
#'
#' @param verbose logical; should information be printed to the console?
#'
#' @import SimDesign stats
#' @return an invisible \code{tibble}/\code{data.frame}-type object of
#' class \code{'Spower'} containing the power results from the
#' simulation experiment
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#'
#' @seealso \code{\link{update}}, \code{\link{SpowerCurve}},
#'   \code{\link{getLastSpower}}, \code{\link{is.CI_within}},
#'   \code{\link{is.outside_CI}}
#'
#' @export
#'
#' @examples
#'
#' ############################
#' # Independent samples t-test
#' ############################
#'
#' # Internally defined p_t.test function
#' args(p_t.test)    # missing arguments required for Spower()
#' # help(p_t.test)  # additional information
#'
#' # p_* functions generate data and return single p-value
#' p_t.test(n=50, d=.5)
#' p_t.test(n=50, d=.5)
#'
#' # test that it works
#' Spower(p_t.test(n = 50, d = .5), replications=10)
#'
#' # also behaves naturally with a pipe
#' p_t.test(n = 50, d = .5) |> Spower(replications=10)
#'
#' \donttest{
#'
#' # Estimate power given fixed inputs (prospective power analysis)
#' out <- Spower(p_t.test(n = 50, d = .5))
#' summary(out)   # extra information
#' as.data.frame(out)  # coerced to data.frame
#'
#' # increase precision (not run)
#' # p_t.test(n = 50, d = .5) |> Spower(replications=30000)
#'
#' # alternatively, increase precision from previous object.
#' #   Here we add 20000 more replications on top of the previous 10000
#' p_t.test(n = 50, d = .5) |>
#'   Spower(replications=20000, lastSpower=out) -> out2
#' out2$REPLICATIONS  # total of 30000 replications for estimate
#'
#' # previous analysis not stored to object, but can be retrieved
#' out <- getLastSpower()
#' out   # as though it were stored from Spower()
#'
#' # Same as above, but executed with multiple cores (not run)
#' p_t.test(n = 50, d = .5) |>
#'    Spower(replications=30000, parallel=TRUE, ncores=2)
#'
#' # Solve N to get .80 power (a priori power analysis)
#' p_t.test(n = NA, d = .5) |>
#'   Spower(power=.8, interval=c(2,500)) -> out
#' summary(out)  # extra information
#' plot(out)
#' plot(out, type = 'history')
#'
#' # total sample size required
#' ceiling(out$n) * 2
#'
#' # same as above, but in parallel with 2 cores
#' out.par <- p_t.test(n = NA, d = .5) |>
#'   Spower(power=.8, interval=c(2,500), parallel=TRUE, ncores=2)
#' summary(out.par)
#'
#' # similar information from pwr package
#' (pwr <- pwr::pwr.t.test(d=.5, power=.80))
#' ceiling(pwr$n) * 2
#'
#' # If greater precision is required and the user has a specific amount of time
#' # they are willing to wait (e.g., 5 minutes) then wait.time can be used. Below
#' # estimates root after searching for 1 minute, and run in parallel
#' #  with 2 cores (not run)
#' p_t.test(n = NA, d = .5) |>
#'   Spower(power=.8, interval=c(2,500), wait.time='1', parallel=TRUE, ncores=2)
#'
#' # Similiar to above for precision improvements, however letting
#' #  the root solver continue searching from an early search history.
#' #  Usually a good idea to increase the maxiter and lower the predCI.tol
#' p_t.test(n = NA, d = .5) |>
#'   Spower(power=.8, interval=c(2,500), lastSpower=out,
#'         maxiter=200, predCI.tol=.008) #starts at last iteration in "out"
#'
#' # Solve d to get .80 power (sensitivity power analysis)
#' p_t.test(n = 50, d = NA) |> Spower(power=.8, interval=c(.1, 2))
#' pwr::pwr.t.test(n=50, power=.80) # compare
#'
#' # Solve alpha that would give power of .80 (criterion power analysis)
#' #    interval not required (set to interval = c(0, 1))
#' p_t.test(n = 50, d = .5) |> Spower(power=.80, sig.level=NA)
#'
#' # Solve beta/alpha ratio to specific error trade-off constant
#' #   (compromise power analysis)
#' out <- p_t.test(n = 50, d = .5) |> Spower(beta_alpha = 2)
#' with(out, (1-power)/sig.level)   # solved ratio
#'
#' # update beta_alpha criteria without re-simulating
#' (out2 <- update(out, beta_alpha=4))
#' with(out2, (1-power)/sig.level)   # solved ratio
#'
#' ##############
#' # Power Curves
#' ##############
#'
#' # SpowerCurve() has similar input, though requires varying argument
#' p_t.test(d=.5) |> SpowerCurve(n=c(30, 60, 90))
#'
#' # solve n given power and plot
#' p_t.test(n=NA, d=.5) |> SpowerCurve(power=c(.2, .5, .8), interval=c(2,500))
#'
#' # multiple varying components
#' p_t.test() |> SpowerCurve(n=c(30,60,90), d=c(.2, .5, .8))
#'
#' ################
#' # Expected Power
#' ################
#'
#' # Expected power computed by including effect size uncertainty.
#' # For instance, belief is that the true d is somewhere around d ~ N(.5, 1/8)
#' dprior <- function(x, mean=.5, sd=1/8) dnorm(x, mean=mean, sd=sd)
#' curve(dprior, -1, 2, main=expression(d %~% N(0.5, 1/8)),
#'       xlab='d', ylab='density')
#'
#' # For Spower, define prior sampler for specific parameter(s)
#' d_prior <- function() rnorm(1, mean=.5, sd=1/8)
#' d_prior(); d_prior(); d_prior()
#'
#' # Replace d constant with d_prior to compute expected power
#' p_t.test(n = 50, d = d_prior()) |> Spower()
#'
#' # A priori power analysis using expected power
#' p_t.test(n = NA, d = d_prior()) |>
#'   Spower(power=.8, interval=c(2,500))
#' pwr::pwr.t.test(d=.5, power=.80) # expected power result higher than fixed d
#'
#'
#' ###############
#' # Customization
#' ###############
#'
#' #   Make edits to the function for customization
#' if(interactive()){
#'     p_my_t.test <- edit(p_t.test)
#'     args(p_my_t.test)
#'     body(p_my_t.test)
#' }
#'
#' # Alternatively, define a custom function (potentially based on the template)
#' p_my_t.test <- function(n, d, var.equal=FALSE, n2_n1=1, df=10){
#'
#'     # Welch power analysis with asymmetric distributions
#'     # group2 as large as group1 by default
#'
#'     # degree of skewness controlled via chi-squared distribution's df
#'     group1 <- rchisq(n, df=df)
#'     group1 <-  (group1 - df) / sqrt(2*df)   # Adjusted mean to 0, sd = 1
#'     group2 <- rnorm(n*n2_n1, mean=d)
#'     dat <- data.frame(group = factor(rep(c('G1', 'G2'),
#'                                      times = c(n, n*n2_n1))),
#'     				  DV = c(group1, group2))
#'     obj <- t.test(DV ~ group, dat, var.equal=var.equal)
#'     p <- obj$p.value
#'     p
#' }
#'
#' # Solve N to get .80 power (a priori power analysis), using defaults
#' p_my_t.test(n = NA, d = .5, n2_n1=2) |>
#'   Spower(power=.8, interval=c(2,500)) -> out
#'
#' # total sample size required
#' with(out, ceiling(n) + ceiling(n * 2))
#'
#' # Solve N to get .80 power (a priori power analysis), assuming
#' #   equal variances, group2 2x as large as group1, large skewness
#' p_my_t.test(n = NA, d=.5, var.equal=TRUE, n2_n1=2, df=3) |>
#'   Spower(power=.8, interval=c(30,100)) -> out2
#'
#' # total sample size required
#' with(out2, ceiling(n) + ceiling(n * 2))
#'
#' # prospective power, can be used to extract the adjacent information
#' p_my_t.test(n = 100, d = .5) |> Spower() -> post
#'
#' ###############################
#' # Using CIs instead of p-values
#' ###############################
#'
#' # CI test returning TRUE if psi0 is outside the 95% CI
#' ci_ind.t.test <- function(n, d, psi0=0, conf.level=.95){
#'   g1 <- rnorm(n)
#'   g2 <- rnorm(n, mean=d)
#'   CI <- t.test(g2, g1, var.equal=TRUE,conf.level=conf.level)$conf.int
#'   is.outside_CI(psi0, CI)
#' }
#'
#' # returns logical
#' ci_ind.t.test(n=100, d=.2)
#' ci_ind.t.test(n=100, d=.2)
#'
#' # simulated prospective power
#' ci_ind.t.test(n=100, d=.2) |> Spower()
#'
#' # compare to pwr package
#' pwr::pwr.t.test(n=100, d=.2)
#'
#' ############################
#' # Equivalence test power using CIs
#' #
#' # H0: population d is outside interval [LB, UB] (not tolerably equivalent)
#' # H1: population d is within interval [LB, UB]  (tolerably equivalent)
#'
#' # CI test returning TRUE if CI is within tolerable equivalence range (tol)
#' ci_equiv.t.test <- function(n, d, tol, conf.level=.95){
#'   g1 <- rnorm(n)
#'   g2 <- rnorm(n, mean=d)
#'   CI <- t.test(g2, g1, var.equal=TRUE,conf.level=conf.level)$conf.int
#'   is.CI_within(CI, tol)
#' }
#'
#' # evaluate if CI is within tolerable interval (tol)
#' ci_equiv.t.test(n=1000, d=.2, tol=c(.1, .3))
#'
#' # simulated prospective power
#' ci_equiv.t.test(n=1000, d=.2, tol=c(.1, .3)) |> Spower()
#'
#' # higher power with larger N (more precision) or wider tol interval
#' ci_equiv.t.test(n=2000, d=.2, tol=c(.1, .3)) |> Spower()
#' ci_equiv.t.test(n=1000, d=.2, tol=c(.1, .5)) |> Spower()
#'
#' ####
#' # superiority test (one-tailed)
#' # H0: population d is less than LB    (not superior)
#' # H1: population d is greater than LB (superior)
#'
#' # set upper bound to Inf as it's not relevant, and reduce conf.level
#' #   to reflect one-tailed test
#' ci_equiv.t.test(n=1000, d=.2, tol=c(.1, Inf), conf.level=.90) |>
#'   Spower()
#'
#' # higher LB means greater requirement for defining superiority (less power)
#' ci_equiv.t.test(n=1000, d=.2, tol=c(.15, Inf), conf.level=.90) |>
#'   Spower()
#'
#' }
Spower <- function(..., power = NA, sig.level=.05, interval,
				   beta_alpha, replications=10000, integer,
				   parallel = FALSE, cl = NULL, packages = NULL,
				   ncores = parallelly::availableCores(omit = 1L),
				   predCI = 0.95, predCI.tol = .01, verbose = TRUE,
				   check.interval = FALSE, maxiter=150, wait.time = NULL,
				   lastSpower = NULL, select = NULL, control = list()){
	if(missing(beta_alpha)) beta_alpha <- NULL
	if(!is.null(cl)) parallel <- TRUE
	control$useAnalyseHandler <- FALSE
	nparent <- control$nparent
	control$nparent <- NULL
	if(is.null(nparent)) nparent <- 1
	pf <- parent.frame(nparent)
	export_funs <- ls(envir = pf)
	if(parallel){
		type <- if(is.null(control$type))
			ifelse(.Platform$OS.type == 'windows', 'PSOCK', 'FORK')
		else control$type
		if(is.null(cl)){
			cl <- parallel::makeCluster(ncores, type=type)
			on.exit(parallel::stopCluster(cl), add = TRUE)
		}
		parallel::clusterExport(cl=cl, export_funs, envir = pf)
		if(verbose)
			message(sprintf("\nNumber of parallel clusters in use: %i", length(cl)))
	}
	packages <- c(packages, 'Spower')
	if(is.na(sig.level)){
		integer <- FALSE
		if(missing(interval)) interval <- c(0, 1)
		stopifnot(interval[2] <= 1 && interval[1] >= 0)
	}
	if(!is.null(wait.time) && maxiter == 150){
		maxiter <- 3000
		predCI.tol <- NULL
	}
	summarise <- Internal_Summarise
	fixed_objects <- list(sig.level=sig.level)
	expr <- match.call(expand.dots = FALSE)$...[[1]]
	expr <- match.call(eval(expr[[1]], envir = pf), expr)
	if(!is.null(expr[-1])){
		pick <- names(which(sapply(expr[-1], \(x){
			ret <- suppressWarnings(try(all(is.na(x)), silent = TRUE))
			if(!is.logical(ret)) ret <- FALSE
			ret
		})))
	} else pick <- character(0)
	fixed_objects$expr <- expr
	fixed_objects$pick <- pick
	fixed_objects$parent_frame <- pf
	stopifnot(is.null(select) || is.character(select))
	fixed_objects$select <- select
	if((is.na(power) + is.na(sig.level) + length(pick)) != 1)
		stop('Exactly *one* argument must be set to \'NA\' in Spower(..., power, sig.level)',
			 call.=FALSE)
	lst_expr <- as.list(expr)[-1]
	if(length(lst_expr))
		lst_expr <- lst_expr[sapply(lst_expr, \(x) is.atomic(x) || is.list(x))]
	conditions <- do.call(SimDesign::createDesign, c(lst_expr, sig.level=sig.level, power=power))
	if(missing(interval)){
		if(is.na(sig.level) || length(fixed_objects$pick))
			stop('Must provide a search interval to solve the missing NA', call.=FALSE)
		interval <- c(NA, NA)
	}
	if(!is.na(power) && !is.na(sig.level) && is.null(beta_alpha)){
		if(missing(integer)){
			integer <- !(has.decimals(interval) || diff(interval) < 5)
			if(!integer && verbose)
				message('\nUsing continuous search interval (integer = FALSE).')
		}
	} else integer <- FALSE
	ret <- if(is.na(power) || !is.null(beta_alpha)){
		conditions$power <- NULL
		if(is.null(beta_alpha))
			summarise <- Internal_Summarise.Full
		seed <- SimDesign::genSeeds(conditions)
		if(!is.null(lastSpower)){
			while(TRUE){
				if(!(seed %in% attr(lastSpower, 'extra_info')$SEED_history)) break
				seed <- SimDesign::genSeeds(conditions)
			}
		}
		tmp <- SimDesign::runSimulation(conditions, replications=replications,
					  analyse=sim_function_aug, summarise=summarise,
					  fixed_objects=fixed_objects, save=FALSE, resume=FALSE,
					  cl=cl, parallel=parallel, ncores=ncores, seed=seed,
					  verbose=verbose, packages=packages, control=control)
		attr(tmp, 'extra_info')$stored_results$sig.level <- NULL
		attr(tmp, 'extra_info')$SEED_history <- seed
		if(!is.null(lastSpower))
			attr(tmp, 'extra_info')$SEED_history <-
			unique(c(seed, attr(lastSpower, 'extra_info')$SEED_history))
		alpha <- 1 - predCI
		pick <- grepl('^power', colnames(tmp))
		pwrnms <- colnames(tmp)[grepl('^power', colnames(tmp))]
		conditions[pwrnms] <- NA
		if(!is.null(lastSpower)){
			attr(tmp, 'extra_info')$stored_results <-
				rbind(attr(tmp, 'extra_info')$stored_results,
					  attr(lastSpower, 'extra_info')$stored_results)
			replications <- tmp$REPLICATIONS + lastSpower$REPLICATIONS
			tmp$REPLICATIONS <- replications
			tmp[pwrnms] <- SimDesign::reSummarise(summarise, results=tmp)[pwrnms]
		}
		CI.lst <- lapply(pwrnms, \(pwrnm){
			CI <- tmp[[pwrnm]] + qnorm(c(alpha/2, predCI+alpha/2)) *
				sqrt((tmp[[pwrnm]] * (1-tmp[[pwrnm]]))/replications)
			CI <- clip_CI(CI)
			CI
		})
		CI <- do.call(rbind, CI.lst)
		rownames(CI) <- pwrnms
		colnames(CI) <- paste0('CI_', c(alpha/2, predCI+alpha/2)*100)
		attr(tmp, 'extra_info')$power.CI <- CI
		attr(tmp, 'extra_info')[c("number_of_conditions", "Design.ID",
								  'save_info', 'functions')] <- NULL
		tmp
	} else {
		SimDesign::SimSolve(conditions, interval=interval,
							analyse=sim_function_aug, save=FALSE,
							summarise=summarise, b=power,
							integer=integer, fixed_objects=fixed_objects,
							cl=cl, parallel=parallel, ncores=ncores,
							verbose=ifelse(verbose, 2, FALSE),
							predCI=predCI, predCI.tol=predCI.tol,
							control=control, check.interval=check.interval,
							maxiter=maxiter, wait.time=wait.time, packages=packages,
							lastSolve=lastSpower)
	}
	if(!is.null(beta_alpha)){
		out <- uniroot(compromise_root, c(.0001, .9999), beta_alpha=beta_alpha,
				sim=ret, Design=conditions, Summarise=Internal_Summarise4Compromise)
		ret$sig.level <- out$root
		ret$power <- 1 - beta_alpha * out$root
		conditions$sig.level <- as.numeric(NA)
		conditions$beta_alpha <- beta_alpha
	}
	attr(ret, 'Spower_extra') <- list(predCI=predCI, conditions=conditions,
							   beta_alpha=beta_alpha, expected=FALSE)
	class(ret) <- c('Spower', class(ret))
	.SpowerEnv$lastSim <- ret
	if(verbose){
		print(ret)
		return(invisible(ret))
	}
	ret
}

sim_function_aug <- function(condition, dat, fixed_objects){
	pick <- fixed_objects$pick
	if(length(pick))
		fixed_objects$expr[pick] <- condition[pick]
	ret <- eval(fixed_objects$expr, envir = fixed_objects$parent_frame)
	if(any(is.logical(ret)))
		ret[is.logical(ret)] <- as.integer(!ret[is.logical(ret)])
	ret
}

#' @rdname Spower
#' @param x object of class \code{'Spower'}
#' @export
print.Spower <- function(x, ...){
	lste <- attr(x, 'Spower_extra')
	cat("\nDesign conditions: \n\n")
	print(lste$conditions)
	if(inherits(x, 'SimSolve')){
		lst <- attr(x, 'roots')[[1]]
		pick <- which(is.na(lste$conditions[1,]))
		cat(sprintf(paste0("\nEstimate of %s: ", if(lst$integer) "%.1f" else "%.3f"),
					names(lste$conditions)[pick],
					x[[pick]]))
		cat(sprintf(paste0("\n%s%% Prediction Interval: ",
						   if(lst$integer) "[%.1f, %.1f]" else "[%.3f, %.3f]", '\n'),
					lste$predCI*100, lst$predCIs_root[1], lst$predCIs_root[2]))
	} else {
		if(!is.null(lste$beta_alpha)){
			cat(sprintf("\nEstimate of Type I error rate (alpha/sig.level): %.3f", x$sig.level))
			alpha <- 1 - lste$predCI
			CI <- x$sig.level + c(qnorm(c(alpha/2, lste$predCI+alpha/2))) *
				sqrt((x$sig.level * (1-x$sig.level))/x$REPLICATIONS)
			CI <- clip_CI(CI)
			cat(sprintf("\n%s%% Confidence Interval: [%.3f, %.3f]\n",
						lste$predCI*100, CI[1], CI[2]))
			power <- x$power
			cat(sprintf("\nEstimate of %spower (1-beta): %.3f",
						if(lste$expected) 'expected ' else "", power))
			CI <- power + c(qnorm(c(alpha/2, lste$predCI + alpha/2))) *
				sqrt((power * (1-power))/x$REPLICATIONS)
			CI <- clip_CI(CI)
			cat(sprintf("\n%s%% Confidence Interval: [%.3f, %.3f]\n",
						lste$predCI*100, CI[1], CI[2]))
		} else {
			CI <- attr(x, 'extra_info')$power.CI
			nms <- if(is.matrix(CI)){
				rownames(CI)
			} else 'power'
			for(nm in nms){
				cat(sprintf("\nEstimate of %s: %.3f", nm, x[[nm]]))
				cat(sprintf("\n%s%% Confidence Interval: [%.3f, %.3f]\n",
							lste$predCI*100, CI[nm,1], CI[nm,2]))
			}
		}
	}
	invisible(NULL)
}


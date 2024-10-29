#' Simulation-based Power Analysis
#'
#' General purpose function that serves as a power-focused wrapper to the
#' \code{SimDesign} package's \code{\link{runSimulation}} and
#' \code{\link{SimSolve}} functions.
#'
#' Five types of power analysis flavors can be performed with \code{Spower},
#' which are triggered based on which supplied inputs as set to missing (\code{NA}):
#'
#' \describe{
#'    \item{A priori}{Solve for a missing sample size component
#'      (e.g., \code{n}) to achieve a specific target power rate
#'      (argument \code{power})}
#'    \item{Post-hoc}{Estimate the power rate given a set of fixed conditions}
#'    \item{Sensitivity}{Solve a missing effect size value as a function of
#'      the other supplied constant components}
#'    \item{Criterion}{Solve the error rate (argument \code{sig.level}) as a
#'      function of the other supplied constant components}
#'    \item{Compromise}{Solve a Type I/Type II error trade-off ratio as a
#'      function of the other supplied constant components and the
#'      target ratio \eqn{q = \beta/\alpha} (argument \code{beta_alpha})}
#' }
#'
#' Post-hoc and compromise analyses utilize the
#' \code{\link[SimDesign]{runSimulation}} function, while the remaining three
#' approaches utilize the stochastic root solving methods in the function
#' \code{\link[SimDesign]{SimSolve}}.
#' See the example below for a demonstration with an independent samples t-test
#' analysis.
#'
#' @param ... a set of conditions to use in the simulation that must match the
#'   arguments in the function \code{sim_function}. Internally these arguments
#'   are passed to either \code{\link{SimSolve}} or
#'  \code{\link{runSimulation}} depending on which element (including
#'  the \code{power} and \code{sig.level} arguments) is set to \code{NA}
#'
#' @param sim function that both creates the data and returns a single
#'   p-value for the analysis of interest
#'
#' @param power power level to use. If set to \code{NA} then the empirical power
#'   will be estimated given the fixed \code{conditions} input
#'   (e.g., for post-hoc power analysis)
#'
#' @param maxiter maximum number of stochastic root-solving iterations
#'
#' @param sig.level alpha level to use. If set to \code{NA} then the empirical
#'   alpha will be estimated given the fixed \code{conditions} input
#'   (e.g., for criterion power analysis)
#'
#' @param interval search interval to use when \code{\link{SimSolve}} is required.
#'   Note that for compromise analyses, where the \code{sig.level} is set to
#'   \code{NA}, if not set explicitly then the interval will default to \code{c(0,1)}
#'
#' @param replications number of replications to use when
#'   \code{\link{runSimulation}} is required
#'
#' @param integer a logical value indicating whether the search iterations
#'   use integers or doubles.
#'   Automatically set to \code{FALSE} if \code{interval} contains
#'   non-integer numbers, as well as when \code{sig.level = NA},
#'   though in general this should be set explicitly
#'
#' @param beta_alpha ratio to use in compromise analyses corresponding to
#'   the Type II errors (beta) over the Type I error (alpha). Ratios greater
#'   than 1 indicate that Type I errors are worse than Type II, while ratios
#'   less than one the opposite. A ratio equal to 1 gives an equal trade-off
#'   between Type I and Type II errors
#'
#' @param parallel for parallel computing for slower simulation experiments
#'   (see \code{\link{runSimulation}} for details)
#'
#' @param cl see \code{\link{runSimulation}}
#'
#' @param ncores see \code{\link{runSimulation}}
#'
#' @param control a list of control parameters to pass to
#'   \code{\link{runSimulation}} or \code{\link{SimSolve}}
#'
#' @param predCI predicting confidence interval level
#'   (see \code{\link{SimSolve}})
#'
#' @param predCI.tol predicting confidence interval consistency tolerance
#'    for stochastic root solver convergence (see \code{\link{SimSolve}}).
#'    Default converges when the power rate CI is consistently
#'    within \code{.01/2} of the target power
#'
#' @param check.interval logical; check the interval range validity
#'   (see \code{\link{SimSolve}})
#'
#' @param verbose logical; should information be printed to the console?
#'
# @param extra_args additional parameters to pass to \code{\link{runSimulation}} or
#   \code{\link{SimSolve}}, specified as a list
#   (e.g., \code{extra_args = list(verbose=FALSE)})
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
#' # Estimate power given fixed inputs (post-hoc power analysis)
#' Spower(p_t.test, n = 50, d = .5)
#'
#' \dontrun{
#'
#' # Same as above, but executed with multiple cores (not run)
#' # Spower(p_t.test, n = 50, d = .5, parallel=TRUE)
#'
#' # Solve N to get .80 power (a priori power analysis)
#' (out <- Spower(p_t.test, n = NA, d = .5, power=.8, interval=c(2,500)))
#' # total sample size required
#' ceiling(out$n) * 2
#'
#' # similar information from pwr package
#' (pwr <- pwr::pwr.t.test(d=.5, power=.80))
#' ceiling(pwr$n) * 2
#'
#' # Solve d to get .80 power (sensitivity power analysis)
#' Spower(p_t.test, n = 50, d = NA, power=.8, interval=c(.1, 2))
#'
#' pwr::pwr.t.test(n=50, power=.80) # compare
#'
#' # Solve alpha that would give power of .80 (criterion power analysis)
#' #    interval not required (set to interval = c(0, 1))
#' Spower(p_t.test, n = 50, d = .5, power=.80, sig.level=NA)
#'
#' # Solve beta/alpha ratio to specific error trade-off constant
#' #   (compromise power analysis)
#' Spower(p_t.test, n = 50, d = .5, beta_alpha = 2)
#'
#'
#' ###############
#' # Customization
#' ###############
#'
#' #   Make edits to the function for customization
#' if(interactive()){
#'     new.p_t.test <- edit(p_t.test)
#'     args(new.p_t.test)
#'     body(new.p_t.test)
#' }
#'
#' # Alternatively, define a custom function (potentially based on the template)
#' new.p_t.test <- function(n, d, var.equal=FALSE, n2_n1=1, df=10){
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
#'     p <- t.test(DV ~ group, dat, var.equal=var.equal)$p.value
#'     p
#' }
#'
#' # Solve N to get .80 power (a priori power analysis), using defaults
#' Spower(p_t.test, n = NA, d = .5, power=.8, interval=c(2,500))
#'
#' # Solve N to get .80 power (a priori power analysis), assuming
#' #   equal variances, group2 2x as large as group1, large skewness
#' (out <- Spower(p_t.test, n = NA, d = .5, var.equal=TRUE, n2_n1=2, df=3,
#'                power=.8, interval=c(2,500)))
#'
#' # total sample size required
#' ceiling(out$n) + ceiling(out$n) * 2
#'
#' # should different alpha level be used given the assumption violations?
#' (TypeI <- Spower(p_t.test, n = 50, d = 0, var.equal=TRUE, n2_n1=2, df=3,
#'                  replications=30000))
#' se <- with(TypeI, sqrt(sig.level * (1-sig.level) / REPLICATIONS))
#' TypeI$sig.level + qnorm(c(.025, .975)) * se  # 95% CI
#' TypeI$power
#'
#'
#' }
Spower <- function(sim, ..., interval, power = NA,
				   sig.level=.05, beta_alpha = NULL,
				   replications=10000, integer,
				   parallel = FALSE, cl = NULL,
				   ncores = parallelly::availableCores(omit = 1L),
				   predCI = 0.95, predCI.tol = .01, verbose = TRUE,
				   check.interval = TRUE, maxiter=150, control = list()){
	fixed_objects <- dots <- list(...)
	dots <- lapply(dots, \(x) if(!is.atomic(x) || length(x) > 1) list(x) else x)
	names(dots) <- names(fixed_objects)
	conditions <- do.call(SimDesign::createDesign, dots)
	stopifnot(nrow(conditions) == 1)
	fixed_objects$ID <- 1
	if(is.na(sig.level) && missing(interval)) interval <- c(0, 1)
	if(is.na(sig.level)) integer <- FALSE
	conditions$sig.level <- fixed_objects$sig.level <- sig.level
	if(!is.na(power)) conditions$power <- power
	if(missing(interval)){
		if(is.na(sig.level) || any(is.na(conditions)))
			stop('Must provide a search interval to solve the missing NA', call.=FALSE)
		interval <- c(NA, NA)
	}
	if(missing(integer)){
		integer <- !has.decimals(interval)
		if(!integer)
			message('Using continuous search interval')
	}
	if(sum(sapply(conditions, \(x) isTRUE(is.na(x))), is.na(power)) != 1)
		stop(c('Exactly one argument for the inputs \'power\', \'sig.level\',',
			   '\n  or the \'...\' list must be set to NA'), call.=FALSE)
	sim_function_aug <- function(condition, dat, fixed_objects){
		pick <- which(sapply(fixed_objects, \(x) all(is.na(x))))
		nm <- names(pick)
		if(length(pick)) fixed_objects[[nm]] <- condition[[nm]]
		do.call(sim,
				fixed_objects[!(names(fixed_objects) %in% c('ID', 'sig.level'))])
	}
	ret <- if(is.na(power) || !is.null(beta_alpha)){
		tmp <- SimDesign::runSimulation(conditions, replications=replications,
					  analyse=sim_function_aug,
					  summarise=Internal_Summarise,
					  fixed_objects=fixed_objects, save=FALSE,
					  cl=cl, parallel=parallel, ncores=ncores,
					  verbose=verbose, control=control)
		tmp <- tmp[,c(names(conditions), 'power')]
		class(tmp) <- c("tbl_df", "tbl", "SimDesign", "data.frame")
		tmp
	} else {
		SimDesign::SimSolve(conditions, interval=interval,
							analyse=sim_function_aug, save=FALSE,
							summarise=Internal_Summarise, b=power,
							integer=integer, fixed_objects=fixed_objects,
							cl=cl, parallel=parallel, ncores=ncores, verbose=verbose,
							predCI=predCI, predCI.tol=predCI.tol,
							control=control, check.interval=check.interval,
							maxiter=maxiter)
	}
	if(!is.null(beta_alpha)){
		out <- uniroot(compromise_root, c(.0001, .9999), beta_alpha=beta_alpha,
				sim=ret, Design=conditions, Summarise=Internal_Summarise4Compromise)
		ret$sig.level <- out$root
		ret$power <- 1 - beta_alpha * out$root
	}
	ret
}

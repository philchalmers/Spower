#' Simulation-based Power Analysis
#'
#' General purpose function that serves as a power-focused wrapper to the
#' \code{SimDesign} package's \code{\link{runSimulation}} and
#' \code{\link{SimSolve}} functions.
#'
#' @param ... a set of conditions to use in the simulation that must match the
#'   arguments in the function \code{sim_function}. Internally these arguments
#'   are passed to either \code{\link{SimSolve}} or
#'  \code{\link{runSimulation}} depending on which element (including
#'  the \code{power} and \code{sig.level} arguments) is set to \code{NA}
#'
#' @param sim_function function that both creates the data and returns a single
#'   p-value for the analysis of interest
#'
#' @param power power level to use. If set to \code{NA} then the empirical power
#'   will be estimated given the fixed \code{conditions} input
#'   (e.g., for post-hoc power analysis)
#'
#' @param sig.level alpha level to use. If set to \code{NA} then the empirical
#'   alpha will be estimated given the fixed \code{conditions} input
#'   (e.g., for criterion power analysis)
#'
#' @param interval search interval to use when \code{\link{SimSolve}} is required
#'
#' @param replications number of replications to use when
#'   \code{\link{runSimulation}} is required
#'
#' @param integer a logical value indicating whether the search iterations
#'   use integers or doubles.
#'   Automatically set to \code{FALSE} if \code{interval} contains
#'   non-integer numbers, though in general this should be set explicitly
#'
#' @param beta_alpha ratio to use in compromise analyses corresponding to
#'   the Type II errors (beta) over the Type I error (alpha). Ratios greater
#'   than 1 indicate that Type I errors are worse than Type II, while ratios
#'   less than one the opposite. A ratio equal to 1 gives an equal trade-off
#'   between Type I and Type II errors
#'
# @param ... additional parameters to pass to \code{\link{runSimulation}} or
#   \code{\link{SimSolve}}
#'
#' @export
#'
#' @examples
#'
#'
#' ############################
#' # Independent samples t-test
#' ############################
#'
#' # Internally defined sim_t.test function
#' args(sim_t.test)    # arguments required for Spower()
#' body(sim_t.test)    # body of the function to generate data + return p-value
#' # help(sim_t.test)  # additional information
#'
#' # Estimate power given fixed inputs (post-hoc power analysis)
#' Spower(N = 50, d = .5, sim_function=sim_t.test)
#'
#' \dontrun{
#'
#' # Solve N to get .80 power (a priori power analysis)
#' out <- Spower(N = NA, d = .5, sim_function=sim_t.test,
#'           	 power=.8, interval=c(2,500))
#' # total sample size required
#' out$N * 2
#'
#' # Solve d to get .80 power (sensitivity power analysis)
#' Spower(N = 50, d = NA, sim_function=sim_t.test,
#' 	   power=.80, interval=c(.1, 2))
#'
#' # Solve alpha that would give power of .80 (criterion power analysis)
#' Spower(N = 50, d = .5, sim_function=sim_t.test,
#' 	   interval=c(.0001, .8), power=.80, sig.level=NA)
#'
#' # Solve beta/alpha ratio to specific trade-off constant
#' #   (compromise power analysis)
#' Spower(N = 50, d = .5, sim_function=sim_t.test, beta_alpha = 2)
#'
#'
#' ###############
#' # Customization
#' ###############
#'
#' #   Make edits to the function for customization
#' if(interactive()){
#'     mysim_t.test <- edit(sim_t.test)
#'     args(mysim_t.test)
#'     body(mysim_t.test)
#' }
#'
#' # Alternatively, define a custom function (potentially based on the template)
#' mysim_t.test <- function(N, d, var.equal=FALSE, group_ratio=1, df=10){
#'
#'     # Welch power analysis with asymmetric distributions
#'     # group1 as large as group2 by default
#'
#'     # degree of skewness controlled via chi-squared distribution's df
#'     group1 <- rchisq(group_ratio*N, df=df)
#'     group1 <-  (group1 - df) / sqrt(2*df)   # Adjusted mean to 0, sd = 1
#'     group2 <- rnorm(N, mean=d)
#'     dat <- data.frame(group = factor(rep(c('G1', 'G2'),
#'                                      times = c(group_ratio*N, N))),
#'     				  DV = c(group1, group2))
#'     p <- t.test(DV ~ group, dat, var.equal=var.equal)$p.value
#'     p
#' }
#'
#' # Solve N to get .80 power (a priori power analysis), using defaults
#' Spower(N = NA, d = .5, sim_function=mysim_t.test,
#'        power=.8, interval=c(2,500))
#'
#' # Solve N to get .80 power (a priori power analysis), assuming
#' #   equal variances, control 2x as large as treatment group, large skewness
#' out <- Spower(N = NA, d = .5, var.equal=TRUE, group_ratio=2, df=3,
#'               sim_function=mysim_t.test, power=.8, interval=c(2,500))
#'
#' # total sample size required
#' out$N * 2 + out$N
#'
#'
#' }
Spower <- function(..., sim_function, interval, power = NA,
				   sig.level=.05, replications=10000, integer,
				   beta_alpha = NULL){
	conditions <- createDesign(...)
	stopifnot(nrow(conditions) == 1)
	conditions$sig.level <- sig.level
	if(missing(interval)){
		if(is.na(sig.level) || any(is.na(conditions)))
			stop('Must provide a search interval to solve the missing NA', call.=FALSE)
		interval <- c(NA, NA)
	}
	if(missing(integer))
		integer <- !has.decimals(interval)
	if(sum(sapply(conditions, is.na), is.na(power)) != 1)
		stop(c('Exactly one argument for the inputs \'power\', \'sig.level\',',
			   '\n  or the \'...\' list must be set to NA'), call.=FALSE)
	fixed_objects <- list()
	sim_function_aug <- function(condition, dat, fixed_objects){
		do.call(sim_function,
				condition[!(names(condition) %in% c('ID', 'sig.level'))])
	}
	ret <- if(is.na(power) || !is.null(beta_alpha)){
		SimDesign::runSimulation(conditions, replications=replications,
					  analyse=sim_function_aug,
					  summarise=Internal_Summarise,
					  fixed_objects=fixed_objects, save=FALSE, )
	} else {
		SimDesign::SimSolve(conditions, interval=interval,
		  		 analyse=sim_function_aug,
		  		 summarise=Internal_Summarise, b=power,
		  		 integer=integer, fixed_objects=fixed_objects, save=FALSE)
	}
	if(!is.null(beta_alpha)){
		out <- uniroot(compromise_root, c(.0001, .9999), beta_alpha=beta_alpha,
				sim=ret, Design=conditions, Summarise=Internal_Summarise4Compromise)
		ret$sig.level <- out$root
		ret$power <- 1 - beta_alpha * out$root
	}
	ret
}

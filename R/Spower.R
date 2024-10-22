#' Simulation-based Power Analysis
#'
#' Description
#'
#' @param conditions a set of conditions to use in the simulation, passed to
#'  either \code{\link{SimSolve}} if exactly one value contains an \code{NA} or to
#'  \code{\link{runSimulation}} if \code{power} is set to \code{NA}
#' @param power power level to use. If set to \code{NA} then the empirical power
#'   will be estimated given the fixed \code{conditions} input
#'   (e.g., for post-hoc power analysis)
#' @param sig.level alpha level to use. If set to \code{NA} then the empirical
#'   alpha will be estimated given the fixed \code{conditions} input
#'   (e.g., for criterion power analysis)
#' @param interval search interval to use when \code{\link{SimSolve}} is required
#' @param sim_function function that both creates the data and returns a single
#'   p-value for the analysis of interest. Function must contain the arguments
#'   \code{condition} to utilize the row in \code{conditions}
#' @param Analyse analyse function (see \code{\link{runSimulation}})
#' @param replications number of replications to use when
#'   \code{\link{runSimulation}} is required
#' @param integer a logical value indicating whether the search iterations
#'   use integers or doubles.
#'   Automatically set to \code{FALSE} if \code{interval} contains
#'   non-integer numbers, though in general this should be set explicitly
#' @param ... additional parameters to pass to \code{\link{runSimulation}} or
#'   \code{\link{SimSolve}}
#'
#' @export
#'
#' @examples
#'
#'
Spower <- function(conditions, sim_function, power, sig.level=.05,
				   interval, replications=10000, integer, compromise.q = NULL,
				   ...){
	conditions$sig.level <- sig.level
	if(missing(interval)) interval <- NA
	if(missing(integer))
		integer <- !has.decimals(interval)
	if(missing(power)) power <- NA
	fixed_objects <- list()
	sim_function_aug <- function(condition, dat, fixed_objects)
		sim_function(condition=condition)
	ret <- if(is.na(power)){
		SimDesign::runSimulation(conditions, replications=replications,
					  analyse=sim_function_aug,
					  summarise=Internal_Summarise,
					  fixed_objects=fixed_objects, ...)
	} else {
		SimDesign::SimSolve(conditions, interval=interval,
		  		 analyse=sim_function_aug,
		  		 summarise=Internal_Summarise, b=power,
		  		 integer=integer, fixed_objects=fixed_objects, ...)
	}
	ret
}

t.test_sim <- function(condition, fixed_objects) {
	Attach(condition)
	group1 <- rnorm(N)
	group2 <- rnorm(N, mean=d)
	dat <- data.frame(group = gl(2, N, labels=c('G1', 'G2')),
					  DV = c(group1, group2))
	p <- t.test(DV ~ group, dat, var.equal=TRUE)$p.value
	p
}

Internal_Summarise <- function(condition, results, fixed_objects) {
	ret <- c(power = EDR(results, alpha = condition$sig.level))
	ret
}

has.decimals <- function(x){
	intx <- as.integer(x)
	any(abs(x - intx) > 0)
}

if(FALSE){

	# estimate power given fixed inputs (post-hoc power analysis)
	Spower(createDesign(N = 50, d=.5),
			   sim_function=t.test_sim)

	# solve N to get 80% power (a priori power analysis)
	Spower(createDesign(N = NA, d = .5),
			   sim_function=t.test_sim,
			   power=.8, interval=c(2,500))

	# solve d to get 80% power (sensitivity power analysis)
	Spower(createDesign(N = 50, d = NA),
			   sim_function=t.test_sim,
			   power=.8, interval=c(.1, 2))

	# solve alpha (criterion power analysis)
	Spower(createDesign(N = 50, d=.5),
			   sim_function=t.test_sim,
			   interval=c(.0001, .8),
			   power=.8, sig.level=NA)

	# beta/alpha given constant q ratio (compromise power analysis)
	Spower(createDesign(N = 50, d=.5),
			   sim_function=t.test_sim,
			   interval=c(.0001, .8),
			   power=.8, sig.level=NA)
}


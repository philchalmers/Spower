#' Run Monte Carlo Simulation Power
#'
#' Description
#'
#' @param conditions a set of conditions to use in the simulation, passed to
#'  either \code{\link{SimSolve}} if exactly one value contains an \code{NA} or to
#'  \code{\link{runSimulation}} if \code{power} is set to \code{NA}
#' @param power power level to use. If set to \code{NA} then
#' @param sig.level alpha level
#' @param interval search interval to use when \code{\link{SimSolve}} is required
#' @param sim_function function that both creates the data and returns a single
#'   p-value for the analysis of interest. Function must contain the arguments
#'   \code{condition} to utilize the row in \code{conditions} and an
#'   argument \code{fixed_objects} for catching any fixed object information
#'   (see \code{\link{runSimulation}} for details)
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
runMCPower <- function(conditions, sim_function, power, sig.level=.05,
					   interval, replications=10000,
					   integer, ...){
	conditions$sig.level <- sig.level
	if(missing(interval)) interval <- NA
	if(missing(integer))
		integer <- !has.decimals(interval)
	if(missing(power)) power <- NA
	sim_function_aug <- function(condition, dat, fixed_objects)
		sim_function(condition=condition, fixed_objects=fixed_objects)
	ret <- if(is.na(power)){
		runSimulation(conditions, replications=replications,
					  analyse=sim_function_aug,
					  summarise=Internal_Summarise, ...)
	} else {
		  SimSolve(conditions, interval=interval,
		  		 analyse=sim_function_aug,
		  		 summarise=Internal_Summarise, b=power,
		  		 integer=integer, ...)
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

	# estimate power given fixed inputs
	runMCPower(createDesign(N = 50, d=.5),
			   sim_function=t.test_sim)

	# solve N to get 80% power
	runMCPower(createDesign(N = NA, d = .5),
			   sim_function=t.test_sim,
			   power=.8, interval=c(2,500))

	# solve d to get 80% power
	runMCPower(createDesign(N = 50, d = NA),
			   sim_function=t.test_sim,
			   power=.8, interval=c(.1, 2))

	# solve alpha
	runMCPower(createDesign(N = 50, d=.5),
			   sim_function=t.test_sim,
			   interval=c(.0001, .8),
			   power=.8, sig.level=NA)
}


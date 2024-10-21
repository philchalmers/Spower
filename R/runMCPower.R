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
#' @param Generate generate function (see \code{\link{runSimulation}})
#' @param Analyse analyse function (see \code{\link{runSimulation}})
#' @param replications number of replications to use when
#'   \code{\link{runSimulation}} is required
#' @param ... additional parameters to pass to \code{\link{runSimulation}} or
#'   \code{\link{SimSolve}}
#'
#' @export
#'
#' @examples
#'
#'
runMCPower <- function(conditions, power, sig.level=.05,
					   interval, Generate, Analyse, replications=10000, ...){
	stopifnot(nrow(conditions) == 1)
	conditions$sig.level <- sig.level
	if(missing(interval)) interval <- NA
	if(missing(power)) power <- NA
	ret <- if(is.na(power)){
		runSimulation(conditions, replications=replications,
					  generate=Generate, analyse=Analyse,
					  summarise=Internal_Summarise, ...)
	} else {
		  SimSolve(conditions, interval=interval,
		  		 generate=Generate, analyse=Analyse,
		  		 summarise=Internal_Summarise, b=power, ...)
	}
	ret
}

t.test_generate <- function(condition, fixed_objects) {
	Attach(condition)
	group1 <- rnorm(N)
	group2 <- rnorm(N, mean=d)
	dat <- data.frame(group = gl(2, N, labels=c('G1', 'G2')),
					  DV = c(group1, group2))
	dat
}

t.test_analyse <- function(condition, dat, fixed_objects) {
	p <- t.test(DV ~ group, dat, var.equal=TRUE)$p.value
	p
}

Internal_Summarise <- function(condition, results, fixed_objects) {
	ret <- c(power = EDR(results, alpha = condition$sig.level))
	ret
}

if(FALSE){

	# solve N to get 80% power
	runMCPower(createDesign(N = NA, d = .5),
			   Generate=t.test_generate,
			   Analyse=t.test_analyse,
			   power=.8, interval=c(2,500))

	# solve d to get 80% power
	runMCPower(createDesign(N = 50, d = NA),
			   Generate=t.test_generate,
			   Analyse=t.test_analyse,
			   power=.8, interval=c(.1, 2), integer=FALSE)

	# estimate power
	runMCPower(createDesign(N = 50, d=.5),
			   Generate=t.test_generate,
			   Analyse=t.test_analyse)

	# solve alpha
	runMCPower(createDesign(N = 50, d=.5), interval=c(.0001, .8),
			   Generate=t.test_generate, power=.8,
			   Analyse=t.test_analyse, sig.level=NA, integer=FALSE)
}


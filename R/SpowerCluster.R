#' Define a parallel cluster object to be used in internal functions
#'
#' This function defines a object that is placed in a relevant internal environment defined in Spower,
#' and stores the object internally so that it can be reused across multiple applications rather than
#' redefining a cluster object for each analysis. Note that this is performed automatically when
#' a call to \code{\link{Spower}} uses \code{parallel = TRUE}. To remove the internal object
#' for the purpose of redefining a cluster use the argument \code{remove = TRUE}.
#'
#' @param spec input that is passed to \code{parallel::makeCluster()}. If no input is given the
#'   maximum number of available local cores minus 1 will be used.
#'   Setting this to NULL will skip a new definition (allows \code{omp_threads} to be used independently)
#' @param ... additional arguments to pass to \code{parallel::makeCluster}
#' @param remove logical; remove previously defined \code{SpowerCluster()}?
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords parallel
#' @export
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'   # use all available cores
#'   SpowerCluster()
#'   SpowerCluster(remove = TRUE)
#'
#' }
#'
#' }
SpowerCluster <- function(spec, remove = FALSE, ...){
	if(!missing(spec) && is.null(spec))
		return(invisible(NULL))
	if(requireNamespace("parallel", quietly = TRUE)){
		if(missing(spec))
			spec <- parallelly::availableCores(omit = 1L)
		if(remove){
			if(is.null(.SpowerClusterEnv$SPOWERCLUSTER)){
				message('There is no visible SpowerCluster() definition.')
				return(invisible(NULL))
			}
			message('Removing previous cluster definition')
			parallel::stopCluster(.SpowerClusterEnv$SPOWERCLUSTER)
			.SpowerClusterEnv$SPOWERCLUSTER <- NULL
			.SpowerClusterEnv$ncores <- 1L
			.SpowerClusterEnv$omp_threads <- 1L
			return(invisible(NULL))
		}
		if(!is.null(.SpowerClusterEnv$SPOWERCLUSTER)){
			return(invisible(NULL))
		}
		message(sprintf('Creating global cluster definition with %i workers ... ',
						ifelse(is.numeric(spec), spec, length(spec))))
		if(is.numeric(spec))
			.SpowerClusterEnv$SPOWERCLUSTER <- parallel::makeCluster(spec, ...)
		else .SpowerClusterEnv$SPOWERCLUSTER <- spec
		.SpowerClusterEnv$ncores <- length(.SpowerClusterEnv$SPOWERCLUSTER)
		parallel::parSapply(cl=.SpowerClusterEnv$SPOWERCLUSTER,
							X=1L:(.SpowerClusterEnv$ncores*2L),
							FUN=function(x) invisible(NULL))
	}
	return(invisible(NULL))
}

.SpowerClusterEnv <- new.env(parent=emptyenv())
.SpowerClusterEnv$ncores <- 1L

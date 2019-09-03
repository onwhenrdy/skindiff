library(foreach)
library(doParallel)
library(Rcpp)

#' @keywords internal
.dskin.sim.helper <- function(parameter, write.to.r, write.to.file)
{
  p <- dskin.toJsonStr(parameter)

  op <- options(digits.secs = 6)
  start.time <- Sys.time()
  res <- .dskin.binding.simulate(p, write.to.r, write.to.file)
  end.time <- Sys.time()
  options(op)

  res <- c(res, list("timestamp" = c("start" = start.time, "end" = end.time),
                     parameter = parameter))
  attr(res, "dskin.class") = "collect.data"
  return (res)
}


#' Run DSkin Simulation
#'
#' Runs a DSkin simulation based on a DSkin Parameter Pack.
#' If a list of configurations is used as input all computations will be run in parallel using the maximum available or needed compute cores
#' of the machine.
#'
#' @param parameter The DSkin Parameter Pack or a list of DSkin Parameter Packs.
#' @param output Output options of the simulation.
#' Options are:
#' \describe{
#'  \item{internal}{Returns the results as an intenal R data structure}
#'  \item{files}{Produces output files in the current working direction}
#' }
#' The maximum number of output options are two and at least one option must be passed.
#'
#' @return DSkin internal R data strucure (base: list) with DSkin class 'collect.data' or named list (file_tag of the configurations as names)
#' of DSkin internal R data structures.
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @useDynLib dskin, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom doParallel stopImplicitCluster
#' @importFrom stats setNames
#'
#' @export
dskin.simulate <- function(parameter, output = c("internal"))
{
  if (length(output) < 1 | length(output) > 2)
  {
    stop("Need at least 1 ouput option and at most 2.")
  }

  write.to.r <- ("internal" %in% output)
  write.to.file <- ("files" %in% output)

  if (!write.to.r & !write.to.file)
  {
    stop("Unknown ouput options defined.")
  }

  is.p.list <- (is.list(parameter) && !.is_dskin_parameter(parameter))

  if (is.p.list)
  {
    if (!all(sapply(parameter, .is_dskin_parameter)))
    {
      stop("At least one parameter in the list is not a DSkin parameter pack.")
    }
  }
  else if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  if (is.p.list)
  {
    names <- sapply(parameter, function(x) {x$log$file_tag})
    registerDoParallel(cores = min(length(parameter), detectCores()))
    i <- 1 # to shut up devtools check
    results <- foreach(i = 1:length(parameter), .final = function(x) setNames(x, names)) %dopar%
    {
      .dskin.sim.helper(parameter[[i]], write.to.r, write.to.file)
    }
    stopImplicitCluster()
    return (results)
  }
  else
  {
    return (.dskin.sim.helper(parameter, write.to.r, write.to.file))
  }
}


#' DSkin Geometry Details
#'
#' Computes geometry details (e.g. minimum step-size, eta) and compartment-baed step-sizes from a DSkin Parameter Pack configuration.
#' If a list of configurations is used as input all computations will be run in parallel using the maximum available or needed compute cores
#' of the machine.
#'
#' @param parameter The DSkin Parameter Pack or a list of DSkin Parameter Packs.
#' @param eta If not Null this will override the mb_eta parameter of (all) DSkin Parameter Pack(s). Valid parameter are 0 < eta < 1.
#' @param resolution If not Null this will override the resolution parameter of (all) DSkin Parameter Pack(s). Valid parameter are resolution > 0.
#'
#' @return DSkin internal R data strucure (base: list) with geometry information for every compartment or named list (file_tag of the configurations as names)
#' of DSkin internal R data structures.
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @useDynLib dskin, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom doParallel stopImplicitCluster
#' @importFrom stats setNames
#'
#' @export
dskin.geometry <- function(parameter, eta = NULL, resolution = NULL)
{
  is.p.list <- (is.list(parameter) && !.is_dskin_parameter(parameter))

  if (is.p.list)
  {
    if (!all(sapply(parameter, .is_dskin_parameter)))
    {
      stop("At least one parameter in the list is not a DSkin parameter pack.")
    }
  }
  else if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  if (!is.null(eta) && (!.isReal(eta) || eta <= 0 | eta >= 1))
  {
    stop("eta must be an integer with 0 < eta < 1.")
  }

  if (!is.null(resolution) && (!.isInteger(resolution) || resolution <= 0))
  {
    stop("resolution must be an integer with resolution > 0.")
  }

  if (is.p.list)
  {
    names <- sapply(parameter, function(x) {x$log$file_tag})
    registerDoParallel(cores = min(length(parameter), detectCores()))
    i <- 1 # to shut up devtools check
    results <- foreach(i = 1:length(parameter), .final = function(x) setNames(x, names)) %dopar%
    {
      param <- parameter[[i]]
      if (!is.null(eta))
      {
        param$sys$mb_eta = eta
      }

      if (!is.null(resolution))
      {
        param$sys$resolution = resolution
      }

      p <- dskin.toJsonStr(param)
      return (.dskin.binding.geometry(p))
    }
    stopImplicitCluster()
    return (results)
  }
  else
  {
    param <- parameter
    if (!is.null(eta))
    {
      param$sys$mb_eta = eta
    }

    if (!is.null(resolution))
    {
      param$sys$resolution = resolution
    }

    p <- dskin.toJsonStr(param)
    return (.dskin.binding.geometry(p))
  }
}

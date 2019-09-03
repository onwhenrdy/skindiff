library(jsonlite)

#' Create a DSkin parameter pack
#'
#' Creates a DSkin parameter pack to configure the DSkin simulation run.
#' This data structure is the fundamental basis for every simulation and post-simulation analysis.
#'
#' @param tag String tag to indicate a project name.
#' @param n.layers The number of barrier layers that shoudl be created.
#' @param layer.names Vector of layer names that should be the same length as n.layers. If NULL, default layer names will be used.
#' @param vehicle.name The name of the vehicle layer.
#' @param sink.name The name of the sink layer.
#' @param use.pk Logical value to indicate whether pharmacokinetic processes shoudl be simulated.
#'
#' @return The DSkin parameter pack.
#'
#' @export
#' @examples
#' p <- dskin.parameter("project_1", 2, c("SC", "DSL"), "TTS", "Blood", TRUE)
dskin.parameter <- function(tag = "Unknown",
                            n.layers = 1,
                            layer.names = NULL,
                            vehicle.name = "Donor",
                            sink.name = "Sink",
                            use.pk = TRUE)
{
  log.p <- .dskin.log_params(tag)
  sys.p <- .dskin.sys_params()

  pk.p <- NULL
  if (use.pk)
  {
    pk.p <- .dskin.pk_params()
  }

  compartment.p <- .dskin.compartment_params(n.layers, layer.names, vehicle.name, sink.name)

  res <- .dskin.compile_params(sys.p, log.p, pk.p, compartment.p)
  attr(res, "dskin.class") = "parameter.pack"
  return (res)
}


#' Parameter Pack to JSON
#'
#' Converts a dskin parameter pack to a JSON string.
#' The string can be read by the DSkin Command Line Tool.
#'
#' @param parameter The DSkin parameter pack (see \code{\link{dskin.parameter}}).
#'
#' @return The JSON representation of the DSkin parameters that can be processed by the DSkin Command Line Tool.
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @importFrom jsonlite toJSON
#' @export
#'
#' @examples
#' p <- dskin.parameter("project_1", 2, c("SC", "DSL"), "Patch", "Blood", TRUE)
#' s <- dskin.toJsonStr(p)
dskin.toJsonStr <-function(parameter)
{
  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  toJSON(parameter, auto_unbox = TRUE, pretty = TRUE)
}


#' Parameter Pack to JSON File
#'
#' Converts a dskin parameter pack to a JSON file.
#'
#' @param parameter The DSkin parameter pack.
#' @param file.name The name of the DSkin JSON file (DSkin configuration parameters) that can be processed by the DSkin Command Line Tool.
#'
#' @return None
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @export
#'
#' @examples
#' p <- dskin.parameter("project_1", 2, c("SC", "DSL"), "Patch", "Blood", TRUE)
#' dskin.toJsonFile(p, "new_cfg.json")
dskin.toJsonFile <-function(parameter, file.name = "dskin_cfg.json")
{
  if (!.isSingleString(file.name))
  {
    stop("file.name must be a single string.")
  }

  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  json.str <- dskin.toJsonStr(parameter)
  cat(json.str, file=file.name)
}


#' DSkin JSON to Parameter Pack
#'
#' Converts a DSkin JSON configuration string into a DSkin pramater pack.
#' \strong{Warning:} The Parameter Pack is not validated.
#'
#' @param json.string The DSkin JSON string.
#'
#' @return The DSkin parameter pack.
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @importFrom jsonlite fromJSON
#' @export
dskin.fromJsonStr <-function(json.string)
{
  res <- fromJSON(json.string)
  attr(res, "dskin.class") = "parameter.pack"

  return(res)
}


#' Compartment Names
#'
#' Returns a list of compartment names from an dskin parameter pack.
#'
#' @param parameter The DSkin parameter pack.
#'
#' @return Vector of compartment names defined in the Parameter Pack including vehicle and sink.
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @export
dskin.compartment.names <- function(parameter)
{
  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  compartment.names <- parameter$compartments$vehicle$name
  layer.list <- parameter$compartments$layers
  if (!is.null(layer.list))
  {
    compartment.names <- c(compartment.names, parameter$compartments$layers$name)
  }
  compartment.names <- c(compartment.names, parameter$compartments$sink$name)

  return (compartment.names)
}


#' List of Compartment Filenames
#'
#' Returns a list of file names associated with mass-over-time profiles defined by a DSkin Parameter Pack.
#'
#' @param parameter The DSkin parameter pack.
#' @param path The path to mass-over-time files (NULL for wd).
#'
#' @return A data frame of file names associated with mass-over-time profiles with named columns (layer names).
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @export
#'
#' @examples
#' p <- dskin.parameter("project_1", 2, c("SC", "DSL"), "Patch", "Blood")
#' file.names <- dskin.mass.filenames(p)
dskin.mass.filenames <-function(parameter, path = NULL)
{
  if (!is.null(path) && !.isSingleString(path))
  {
    stop("path must be a single string.")
  }

  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  file.tag <- parameter$log$file_tag
  mass.tag <- parameter$log$mass_file_postfix
  is.zipped <- parameter$log$mass_file_gzip

  c.names <- dskin.compartment.names(parameter)
  # first is vehicle, last is sink
  # check for logged compartments
  has.log <- c(parameter$compartments$vehicle$log,
               parameter$compartments$layers$log,
               parameter$compartments$sink$log)
  c.names <- c.names[has.log]

  c.files <- sapply(c.names, function(x)
  {
    tmp <- paste(path, file.tag, "_" , x, "_", mass.tag, ".dat", sep="")
    if (is.zipped)
    {
      tmp <- paste(tmp, ".gz",sep="")
    }
    return (tmp)
  }
  )

  return (c.files)
}


#' List of Compartment Filenames
#'
#' Returns a list of file names associated with concentration-depth profiles (CDP) defined by a DSkin Parameter Pack.
#'
#' @param parameter The DSkin parameter pack.
#' @param path The path to CDP files (NULL for wd).
#'
#' @return A data frame of file names associated with CDPs with named columns (layer names).
#' @importFrom utils head
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @export
#'
#' @examples
#' p <- dskin.parameter("project_1", 2, c("SC", "DSL"), "Patch", "Blood")
#' file.names <- dskin.conc.filenames(p)
dskin.conc.filenames <-function(parameter, path = NULL)
{
  if (!is.null(path) && !.isSingleString(path))
  {
    stop("path must be a single string.")
  }

  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  file.tag <- parameter$log$file_tag
  conc.tag <- parameter$log$cdp_file_postfix
  is.zipped <- parameter$log$cdp_file_gzip

  c.names <- dskin.compartment.names(parameter)
  # first is vehicle, last is sink
  # check for logged compartments
  has.log <- c(parameter$compartments$vehicle$log_cdp,
               parameter$compartments$layers$log_cdp,
               FALSE)
  c.names <- c.names[has.log]

  c.files <- sapply(c.names, function(x)
  {
    tmp <- paste(path, file.tag, "_" , x, "_", conc.tag, ".dat", sep="")
    if (is.zipped)
    {
      tmp <- paste(tmp, ".gz",sep="")
    }
    return (tmp)
  }
  )

  return (c.files)
}


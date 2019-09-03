########## helper functions


#' @keywords internal
.is_dskin_parameter <- function(ds)
{
  if (is.null(ds))
  {
    return (FALSE)
  }
  test <- attr(ds, "dskin.class")
  if (is.null(test))
  {
    return (FALSE)
  }

  return (test == "parameter.pack")
}


#' @keywords internal
.is_dskin_ss_fit_data <- function(ds)
{
  if (is.null(ds))
  {
    return (FALSE)
  }
  test <- attr(ds, "dskin.class")
  if (is.null(test))
  {
    return (FALSE)
  }

  return (test == "ss.fit.data")
}


#' @keywords internal
.is_dskin_mass_data <- function(ds)
{
  if (is.null(ds))
  {
    return (FALSE)
  }

  test <- attr(ds, "dskin.class")
  if (is.null(test))
  {
    return (FALSE)
  }

  return (test == "mass.data")
}


#' @keywords internal
.is_dskin_ds <- function(ds)
{
  if (is.null(ds))
  {
    return (FALSE)
  }

  test <- attr(ds, "dskin.class")
  if (is.null(test))
  {
    return (FALSE)
  }

  return (test == "collect.data")
}


#' @keywords internal
.is_dskin_conc_data <- function(ds)
{
  if (is.null(ds))
  {
    return (FALSE)
  }

  test <- attr(ds, "dskin.class")
  if (is.null(test))
  {
    return (FALSE)
  }

  return (test == "conc.data")
}

#' @keywords internal
.check_for_files_access <- function(files)
{
  res <- list(access = TRUE, msg = "")

  f.exists <- file.exists(files)
  if (!all(f.exists))
  {
    res[[1]] = FALSE
    res[[2]] = paste("Cannot find files:", paste(files[!f.exists], collapse = ", "))
  }
  else
  {
    # read permissions
    f.access <- file.access(files, 4)
    if (!all(f.access))
    {
      res[[1]] = FALSE
      res[[2]] = paste("Cannot access files:", paste(files[!f.access], collapse = ", "))
    }
  }

  return(res)
}


#' @keywords internal
.isSingleString <- function(input)
{
  if (is.null(input))
  {
    return (FALSE)
  }

  return (is.character(input) && length(input) == 1)
}

#' @keywords internal
.isInteger <- function(input)
{
  if (length(input) != 1)
  {
    return (FALSE)
  }

  if (is.null(input))
  {
    return (FALSE)
  }

  as.int <- as.integer(input)
  if (is.na(as.int))
  {
    return (FALSE)
  }


  all(input == as.integer(input))
}

#' @keywords internal
.isBool <- function(input)
{
  if (length(input) != 1)
  {
    return (FALSE)
  }

  if(is.na(input))
  {
    return (FALSE)
  }

  return(is.logical(input))
}

#' @keywords internal
.isReal <- function(input)
{
  if (length(input) != 1)
  {
    return (FALSE)
  }

  if (is.null(input))
  {
    return (FALSE)
  }

  as.int <- as.integer(input)
  if (is.na(as.int))
  {
    return (FALSE)
  }

  all(input == as.numeric(input))
}


#' @keywords internal
#' @importFrom stats lm
.ss_finder_fit <-function(times, masses, max.window.size = 1, diff.lag = 1, threshold = 0.9)
{
  t.size <- length(times)
  if (t.size < (max.window.size * 2 +1))
  {
    warning("Not enough time points for steady-state calculations found.")
    return(NULL)
  }

  slope <- NA
  slope.time <- NA
  s.idx <- NA

  tmp <- diff(masses, lag= diff.lag)
  max.slope.idx = which(tmp == max(tmp))[1]
  if (!is.na(max.slope.idx) && !is.nan(max.slope.idx))
  {
    # construct fitting window
    min.idx <- max(max.slope.idx - max.window.size, 1)
    max.idx <- min(max.slope.idx + max.window.size, t.size)


    t_s <- times[min.idx:max.idx]
    m_s <- masses[min.idx:max.idx]
    options(warn=-1)
    res <- lm(m_s ~ t_s)
    result <- summary(res)
    options(warn=0)
    if (!is.nan(result$r.squared) && result$r.squared >= threshold)
    {
      slope <- result$coefficients[2]
      s.idx <- ceiling((max.idx - min.idx)/2) + min.idx
      slope.time <- times[s.idx]
    }
  }

  if (is.na(slope))
  {
    warning("Could not find steady-state.")
  }

  ret <- list (slope = slope,
               time = slope.time,
               slope.idx = s.idx)
  return (ret)
}


#' Convert masses
#'
#' Converts one or more values of a certain mass unit into another mass unit.
#'
#' @param val The value [int or real].
#' @param val_unit The unit of the value [mg/ug/ng].
#' @param target_unit The unit the value should be converted into [mg/ug/ng].
#'
#' @return Converted masses (single value or compatible data structure)
#' @export
#'
#' @examples
#' # convert 10 ng to mg
#' dskin.convert.masses(10, "ng", "mg")
dskin.convert.masses <-function(val, val_unit, target_unit)
{
  cov.list = list(mg = 1E6, ug = 1E3, ng = 1)

  # convert first to ng
  val.conv = cov.list[[val_unit]]
  if (is.null(val.conv))
  {
    stop ("Unknown val_unit")
  }

  base.val <-  val.conv * val

  tar.conv = cov.list[[target_unit]]
  if (is.null(tar.conv))
  {
    stop ("Unknown target_unit")
  }

  return (base.val / tar.conv)
}


#' Convert times
#'
#' Converts one or more values of a certain time unit into another time unit.
#'
#' @param val The value [int or real].
#' @param val_unit The unit of the value [sec/min/h].
#' @param target_unit The unit the value should be converted into [sec/min/h].
#'
#' @return Converted times (single value or compatible data structure)
#' @export
#'
#' @examples
#' # convert 10 sec to min
#' dskin.convert.times(10, "sec", "min")
dskin.convert.times <-function(val, val_unit, target_unit)
{
  cov.list = list(h = 3600, min = 60, sec = 1)

  # convert first to sec
  base.val <- val
  val.conv = cov.list[[val_unit]]
  if (is.null(val.conv))
  {
    stop ("Unknown val_unit")
  }

  base.val <-  val.conv * val

  tar.conv = cov.list[[target_unit]]
  if (is.null(tar.conv))
  {
    stop ("Unknown target_unit")
  }

  return (base.val / tar.conv)
}


#' Convert spaces
#'
#' Converts one or more values of a certain space unit into another space unit.
#'
#' @param val The value [int or real].
#' @param val_unit The unit of the value [um/cm/dm].
#' @param target_unit The unit the value should be converted into [um/cm/dm].
#'
#' @return Converted spaces (single value or compatible data structure)
#' @export
#'
#' @examples
#' # convert 10 dm to um
#' dskin.convert.spaces(10, "dm", "um")
dskin.convert.spaces <-function(val, val_unit, target_unit)
{
  cov.list = list(dm = 1E5, cm = 1E4, um = 1)

  # convert first to um
  base.val <- val
  val.conv = cov.list[[val_unit]]
  if (is.null(val.conv))
  {
    stop ("Unknown val_unit")
  }

  base.val <-  val.conv * val

  tar.conv = cov.list[[target_unit]]
  if (is.null(tar.conv))
  {
    stop ("Unknown target_unit")
  }

  return (base.val / tar.conv)
}



#' @keywords internal
.dskin.check_in_path <-function(filename)
{
  if (!.isSingleString(filename))
  {
    stop("Function needs a single single string as input.")
  }

  path <- Sys.which(filename)[[1]]
  return (path != "" & !is.null(path))
}


#' @keywords internal
.temp_unzip <- function(filename, fun, ...){
  BFR.SIZE <- 1e7
  if (!file.exists(filename)) {
    stop("No such file: ", filename);
  }
  if (!is.function(fun)) {
    stop(sprintf("Argument 'fun' is not a function: %s", mode(fun)));
  }
  temp_dir <- tempdir()
  # test if it's zip
  files_in_zip <- try(utils::unzip(filename, list = TRUE)$Name, silent = TRUE)
  if (class(files_in_zip) == "character") {
    # hidden files can be ignored: starting with ., ending with $, __MACOSX folder
    visible_files <- files_in_zip[!grepl("((^__MACOSX\\/.*)|(^\\..*)|(^.*\\$$))",
                                         files_in_zip)]
    # will not continue for multiple non-hidden files since behavior is not well defined.
    if(length(visible_files)>1) {
      stop(paste0("Zip file contains multiple visible files:\n",
                  paste0("    ", visible_files, collapse = "\n")))
    }
    if(length(visible_files) == 0) { stop("\n  No visible file found in Zip file")}
    # proceed with single non-hidden file
    utils::unzip(filename, files = visible_files[1], exdir = temp_dir, overwrite = TRUE)
    dest_file <- file.path(temp_dir, visible_files[1])
  } else {
    dest_file <- tempfile()
    # Setup input and output connections
    inn <- gzfile(filename, open = "rb")
    out <- file(description = dest_file, open = "wb")
    # Process
    nbytes <- 0
    repeat {
      bfr <- readBin(inn, what=raw(0L), size=1L, n=BFR.SIZE)
      n <- length(bfr)
      if (n == 0L) break;
      nbytes <- nbytes + n
      writeBin(bfr, con=out, size=1L)
      bfr <- NULL  # Not needed anymore
    }
    close(inn)
    close(out)
  }
  # call fun with temp file
  res <- fun(dest_file, ...)
  file.remove(dest_file)
  return(res)
}


#' Calulcate Initial Concentration
#'
#' Calculates an initial concentration from compartment dimensions and loading dose assuming a homogenous
#' distribution of mass inside the compartment.
#'
#' @param dose.in.mg The loading dose in mg.
#' @param height.um  The compartment height in micrometer.
#' @param A.sqcm The vehicle application area in \eqn{cm^{2}}{cm^2}
#'
#' @return The homogeneous concentration inside the compartment in mg/ml.
#' @export
#'
#' @examples
#' # Initial concentration of 10 mg inside a compartment of 100 micrometer x 15 cm^2
#' dskin.c.init(10, 100, 15)
dskin.c.init <-function(dose.in.mg, height.um, A.sqcm)
{
  h.cm <- height.um * 1E-4
  return (dose.in.mg / (A.sqcm * h.cm))
}


#' DSkin Runtime String
#'
#' Produces a human readable string of the runtime of a DSkin simulation
#'
#' @param result The DSkin result data structure (dskin class = 'collect.data').
#'
#' @return A human readable string of the runtime.
#'
#' @seealso \code{\link{dskin.simulate}}
#'
#' @export
dskin.runtime.str <-function(result)
{
  if (!.is_dskin_ds(result))
  {
    stop("result is no DSkin data structure")
  }

  format(result$timestamp[2] - result$timestamp[1])
}


#' Parameter Pack Log Settings
#'
#' Shortcut function to set log settings of a DSkin Parameter Pack.
#'
#' @param parameter The input DSkin Parameter Pack.
#' @param mass.logging Logical value if mass-over-time logging should be enabled (TRUE) or disabled (FALSE) for all compartments.
#' @param cdp.logging Logical value if concentration-depth profile logging should be enabled (TRUE) or disabled (FALSE) for all compartments.
#' @param mass.interval Logging inverval in simulated minutes (>0) for mass-over-time logging.
#' @param cdp.interval Logging inverval in simulated minutes (>0) for concentration-depth profile logging.
#'
#' @return The modified DSkin Parameter Pack.
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @export
dskin.log.settings <-function(parameter, mass.logging = TRUE, cdp.logging = TRUE,
                            mass.interval = 1, cdp.interval = 1)
{
  if(!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin Parameter Pack.")
  }

  if (!.isInteger(mass.interval) || mass.interval < 1)
  {
    stop("mass.interval must be an integer with mass.interval > 0")
  }

  if (!.isInteger(cdp.interval) || cdp.interval < 1)
  {
    stop("cdp.interval must be an integer with mass.interval > 0")
  }

  if (!.isBool(mass.logging))
  {
    stop("mass.logging must be logical.")
  }

  if (!.isBool(cdp.logging))
  {
    stop("cdp.logging must be logical.")
  }

  p <- parameter
  p$log$mass_log_interval = mass.interval
  p$log$cdp_log_interval = cdp.interval

  p$compartments$vehicle$log = mass.logging
  p$compartments$vehicle$log_cdp = cdp.logging

  p$compartments$sink$log = mass.logging

  p$compartments$layers$log = mass.logging
  p$compartments$layers$log_cdp = cdp.logging

  return (p)
}

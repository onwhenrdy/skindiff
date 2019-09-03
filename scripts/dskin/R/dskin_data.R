library(DescTools)
library(data.table)

#' Read Mass-over-time Data
#'
#' Reads mass-over-time data from files produced by the DSkin Commandline Tool.
#'
#' @param parameter The DSkin parameter pack (see \code{\link{dskin.parameter}}).
#' @param path The path to the mass-over-time files (NULL for wd).
#'
#' @return Data frame with named compartments and DSkin class 'mass.data'.
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @importFrom utils read.csv
#' @export
dskin.read.massData <- function(parameter, path = NULL)
{
  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  c.files <- dskin.mass.filenames(parameter, path)
  check.access <- .check_for_files_access(c.files)
  if (!check.access[[1]])
  {
    stop(check.access[[2]])
  }

  is.zipped <- parameter$log$mass_file_gzip

  unzip.fn <- function(x) {x}
  if (is.zipped)
  {
    unzip.fn = gzfile
  }

  tmp <- sapply(c.files, function(x) {read.csv(x, header = TRUE, sep="\t")})
  tmp <- rbind(tmp, file.name = c.files)
  tmp <- as.data.frame(tmp)
  attr(tmp, "dskin.class") = "mass.data"
  return (tmp)
}


#' Read Concentration-depth Profiles
#'
#' Reads concentration-depth profiles from files produced by the DSkin Commandline Tool.
#'
#' @param parameter The DSkin parameter pack (see \code{\link{dskin.parameter}}).
#' @param path The path to the concentration-depth profile files (NULL for wd).
#'
#' @return A list with named compartments and DSkin class 'conc.data'.
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom data.table fread
#' @export
dskin.read.concData <- function(parameter, path = NULL)
{
  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a dskin parameter pack.")
  }

  c.files <- dskin.conc.filenames(parameter, path)
  check.access <- .check_for_files_access(c.files)
  if (!check.access[[1]])
  {
    stop(check.access[[2]])
  }

  is.zipped <- parameter$log$cdp_file_gzip

  read.fn <- function(x, ...) {fread(x, ...)}
  if (is.zipped)
  {
    read.fn <- function(x, ...) {.temp_unzip(x, fread, ...)}
  }

  tmp <- lapply(c.files, function(x) {read.fn(x, sep="\t", header = FALSE, data.table = FALSE)})
  tmp <- rbind(data = tmp, file.name = c.files)
  tmp <- as.data.frame(tmp)
  tmp <- as.list(tmp)
  for (comp in names(tmp))
  {
    # add list elements for conc-data, filename, times, and space
    tmp[[comp]] = c(tmp[[comp]],
                    tail(tmp[[comp]]$data[1], -1),
                    as.data.frame(tail(t(tmp[[comp]]$data[1,]), -1)))
    names(tmp[[comp]]) <- c("data", "file.name", "time", "x")

    # remove first row and first column from data (space and time)
    tmp[[comp]]$data <- tmp[[comp]]$data[-c(1),]
    tmp[[comp]]$data[, 1] <- NULL

    # transpose ->
    # columns are times
    # rows are space
    tmp[[comp]]$data <- as.data.frame(t(tmp[[comp]]$data))
    colnames(tmp[[comp]]$data) <- as.character(tmp[[comp]]$time)
    rownames(tmp[[comp]]$data) <- as.character(tmp[[comp]]$x)
  }
  attr(tmp, "dskin.class") = "conc.data"

  return (tmp)
}


#' Compartmental Mass Analysis
#'
#' Calculates maximum mass, time to maximum mass, minimum mass and time to minimum mass for all logged compartments.
#'
#' @param parameter The DSkin Parameter Pack.
#' @param r.result DSkin result data structure ('collect.data') by the R-version of DSkin. If NULL, data is read from produced files
#' defined by the Parameter Pack.
#' @param path The path to the mass-over-time files (NULL for wd).
#' @param time.unit The unit of time used for the overview [sec,min,h].
#' @param mass.unit The mass unit for all layer compartments [ng,ug,mg].
#'
#' @return Data.frame with compartmental names as column names and metrics as row names.
#' The following attributes are attached to the data.frame:
#' \describe{
#'  \item{mass.unit}{The unit used for masses}
#'  \item{time.unit}{The unit used for times}
#' }
#'
#' @seealso \code{\link{dskin.parameter}}
#' @seealso \code{\link{dskin.simulate}}
#'
#' @export
dskin.mass.analyze <- function(parameter,
                               r.result = NULL,
                               path = NULL,
                               time.unit = "h",
                               mass.unit = "mg")
{
  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  # read data
  if (is.null(r.result))
  {
    mass.data <- dskin.read.massData(parameter, path)
  }
  else
  {
    if (!.is_dskin_ds(r.result))
    {
      stop("r.result is not a DSkin data structure.")
    }

    mass.data <- r.result$masses
  }

  m.max = c()
  t.max = c()
  m.min = c()
  t.min = c()

  sink.name = parameter$compartments$sink$name
  for (comp in names(mass.data))
  {
    dat = mass.data[[comp]]

    # max/min masses/times
    m.max.local = max(dat$mass)
    t.max.local = dat$time[which.max(dat$mass)]
    m.min.local = min(dat$mass)
    t.min.local = dat$time[which.min(dat$mass)]

    if (comp == sink.name)
    {
      V = parameter$compartments$sink$Vd
      m.max.local <- m.max.local * V
      m.min.local <- m.min.local * V

    }

    m.max = c(m.max, m.max.local)
    t.max = c(t.max, t.max.local)
    m.min = c(m.min, m.min.local)
    t.min = c(t.min, t.min.local)
  }

  unit <- parameter$log$scaling
  # all units are in mg / cm^2 / cm / hours or combination of all
  df <- data.frame(rbind(max.mass = dskin.convert.masses(m.max, unit, mass.unit),
                         max.time = dskin.convert.times(t.max, "min", time.unit),
                         min.mass = dskin.convert.masses(m.min, unit, mass.unit),
                         min.time = dskin.convert.times(t.min, "min", time.unit)))
  names(df) <- names(mass.data)
  attr(df, "mass.unit") = mass.unit
  attr(df, "time.unit") = time.unit
  return (df)
}


#' Calculate Compartment AUC
#'
#' Calculates the Area under the Curve (AUC) in units of mass.unit * time.unit / ml for every logged compartment.
#' AUC.last calculation is only performed if pharamacokinietics is enabled in the sink compartment and
#' assumes no transport from skin layers to the sink compartment at the last simulated time point.
#'
#' @param parameter The DSkin Parameter Pack.
#' @param r.result DSkin result data structure ('collect.data') by the R-version of DSkin. If NULL, data is read from produced files
#' defined by the Parameter Pack.
#' @param path The path to the mass-over-time files (NULL for wd).
#' @param time.unit The unit of time used for the overview [sec,min,h].
#' @param mass.unit The mass unit for all layer compartments [ng,ug,mg].
#' @param scale.with.cross.section Logical value to indicate if layer cross-section is taken into account for
#' volume calculation of the different skin layers.
#'
#' @return Data.frame with compartmental names as column names and metrics as row names.
#' The following attributes are attached to the data.frame:
#' \describe{
#'  \item{mass.unit}{The unit used for masses}
#'  \item{time.unit}{The unit used for times}
#'  \item{vol.unit}{The unit used for volumes. Currently always ml, but might change for future releases.}
#' }
#'
#' @seealso \code{\link{dskin.parameter}}
#' @seealso \code{\link{dskin.simulate}}
#'
#' @importFrom DescTools AUC
#'
#' @export
dskin.AUC.analyze <- function(parameter,
                              r.result = NULL,
                              path = NULL,
                              time.unit = "h",
                              mass.unit = "mg",
                              scale.with.cross.section = FALSE)
{
  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  # read data
  if (is.null(r.result))
  {
    mass.data <- dskin.read.massData(parameter, path)
  }
  else
  {
    if (!.is_dskin_ds(r.result))
    {
      stop("r.result is not a DSkin data structure.")
    }

    mass.data <- r.result$masses
  }

  AUC.data <- rep(NA, ncol(mass.data))
  AUC.last.data <- rep(NA, ncol(mass.data))

  # in cm^2
  app.area = parameter$compartments$vehicle$app_area
  v.name <- parameter$compartments$vehicle$name
  s.name <- parameter$compartments$sink$name
  dat.unit <- parameter$log$scaling
  for (comp in names(mass.data))
  {
    idx <- which(names(mass.data) == comp)
    dat = mass.data[[comp]]
    auc <- AUC(dat$time, dat$mass)
    auc <- dskin.convert.masses(auc, dat.unit, mass.unit)
    auc <- dskin.convert.times(auc, "min", time.unit)

    if (comp == s.name)
    {
      # sink data is always in mass / ml

      if (parameter$PK$enabled)
      {
        # AUC last
        # t_half is in hours
        t.half <- dskin.convert.times(parameter$PK$t_half, "h", time.unit)
        k.el <- log(2)/t.half
        auc.last <- dskin.convert.masses(tail(dat$mass, 1), dat.unit, mass.unit) / k.el
        AUC.last.data[idx] = auc.last
      }
    }
    else if (comp == v.name)
    {
      # in cm^3 = ml
      v.vol = app.area * parameter$compartments$vehicle$h * 1E-4
      auc = auc / v.vol
    }
    else
    {
      r.idx <- which(parameter$compartments$layers$name == comp)
      h <- parameter$compartments$layers$h[r.idx]
      cs <- ifelse(scale.with.cross.section, parameter$compartments$layers$cross_section[r.idx], 1)
      # in cm^3 = ml
      v.vol = app.area * h * 1E-4 * cs
      auc = auc / v.vol
    }

    AUC.data[idx] = auc
  }

  df <- data.frame(rbind(AUC = AUC.data,
                         AUC.last = AUC.last.data))
  names(df) <- names(mass.data)
  attr(df, "mass.unit") = mass.unit
  attr(df, "time.unit") = time.unit
  attr(df, "vol.unit") = "ml"
  return (df)
}


#' Steady-state Analyzer
#'
#' Function to find and analyze the steady-state of a simulated curve in the sink compartment.
#' Mass-over-time data of the sink compartment is needed to perform the caluclations.
#'
#' @param parameter The DSkin Parameter Pack.
#' @param r.result DSkin result data structure ('collect.data') by the R-version of DSkin. If NULL data is read from files defined by the Parameter Pack.
#' @param path The path to CDP files (NULL for wd).
#' @param mass.unit The mass unit for all layer compartments [ng,ug,mg].
#' @param time.unit The unit of time used for the overview [sec,min,h].
#' @param space.unit The unit of space used of the overview [um, cm, dm]
#'
#' @return Data frame with file_tag as column name and the following rows:
#' \describe{
#'  \item{j.ss}{Steady-state flux in mass.unit / (time.unit * space.unit**2) (e.g. mg/h*cm^2)}
#'  \item{time.ss}{Time until steady-state was reached (might be too long due to small variations in accumulated mass; e.g. in h)}
#'  \item{kp}{Permeability coefficient in space.unit / time.unit (e.g. cm/h)}
#'  \item{time.lag}{Lag-time in time.unit (e.g. h)}
#' }
#'
#' The following attributes are attached:
#' \describe{
#'  \item{mass.unit}{The unit used for masses}
#'  \item{time.unit}{The unit used for times}
#'  \item{space.unit}{The unit used for space units.}
#' }
#'
#' @seealso \code{\link{dskin.parameter}}
#' @seealso \code{\link{dskin.simulate}}
#'
#' @export
dskin.ss.analyze <- function(parameter,
                             r.result = NULL,
                             path = NULL,
                             mass.unit = "mg",
                             time.unit = "h",
                             space.unit = "cm")
{
  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  # read data
  if (is.null(r.result))
  {
    mass.data <- dskin.read.massData(parameter, path)
  }
  else
  {
    if (!.is_dskin_ds(r.result))
    {
      stop("r.result is not a DSkin data structure.")
    }

    mass.data <- r.result$masses
  }

  sink.name = parameter$compartments$sink$name
  if (!(sink.name %in% names(mass.data)))
  {
    stop("Could not find sink compartment in data. Did you forgot to log the compartment?")
  }

  dat = mass.data[[sink.name]]
  # mass unit
  unit <- parameter$log$scaling
  # A is in cm^2
  A = parameter$compartments$vehicle$app_area
  # V is in ml
  V = parameter$compartments$sink$Vd
  # in mg/ml
  c.init <- parameter$compartments$vehicle$c_init

  j.ss  = NA
  t.ss  = NA
  kp    = NA
  t.lag = NA
  ss.find <- .ss_finder_fit(dat$time, dat$mass)
  if (!is.na(ss.find$slope))
  {
    # slope unit: unit/(ml * min)
    # J_ss = slope/ (V * A)
    # J_ss is in unit / (min * cm^2)
    j.ss <- ss.find$slope / (V * A)
    # target is mass.unit / (time.unit * space.unit^2)
    j.ss <- dskin.convert.masses(j.ss, unit, mass.unit)
    j.ss <- j.ss / dskin.convert.times(1, "min", time.unit)
    j.ss <- j.ss / dskin.convert.times(1, "min", time.unit)
    j.ss <- j.ss / (dskin.convert.spaces(1, "cm", space.unit)**2)

    # kp is in space.unit/time.unit
    c.init <- dskin.convert.masses(c.init, "mg", mass.unit)
    c.init <- c.init / (dskin.convert.spaces(1, "cm", space.unit))**3
    kp <- j.ss / c.init

    # in time.unit
    t.ss <- dskin.convert.times(ss.find$time, "min", time.unit)
    t.lag <- dskin.convert.times( - dat$mass[ss.find$slope.idx]/ss.find$slope + ss.find$time, "min", time.unit)
  }

  df <- data.frame(rbind(j.ss,
                         time.ss = t.ss,
                         kp,
                         time.lag = t.lag))
  names(df) <- parameter$log$file_tag
  attr(df, "mass.unit") <- mass.unit
  attr(df, "space.unit") <- space.unit
  attr(df, "time.unit") <- time.unit
  attr(df, "dskin.class") <- "ss.fit.data"
  return (df)
}






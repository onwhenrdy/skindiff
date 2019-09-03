library(foreach)
library(doParallel)

#' Overview Mass-over-time Plots
#'
#' Produces an overview plot for all mass compartments and an optional overall system mass-over-time summary.
#'
#' @param parameter The DSkin Parameter Pack.
#' @param r.result DSkin result data structure ('collect.data') by the R-version of DSkin. If NULL, data is read from produced files
#' defined by the Parameter Pack.
#' @param path The path to the mass-over-time files (NULL for wd).
#' @param max.cols The maximum number of columns in the plot.
#' @param plot.total Logical value to define if total system mass-over-time should be plotted.
#' @param total.main String of the system (plot.total) mass-over-time plot heading.
#' @param time.unit The unit of time used for the plot [sec,min,h].
#' @param layer.unit The mass unit for all layer compartments [ng,ug,mg].
#' @param sink.unit The mass unit for the sink compartment [ng,ug,mg].
#' @param mass.range The range of plottable mass (NULL for automatic min/max).
#' @param time.range The time range for the plots (NULL for automatic min/max).
#' @param ... Additional parameter for plot directive.
#'
#' @return None
#'
#' @seealso \code{\link{dskin.parameter}}
#' @seealso \code{\link{dskin.simulate}}
#'
#' @importFrom graphics plot
#' @importFrom graphics par
#' @importFrom graphics par
#' @export
dskin.plot.masses <- function(parameter,
                              r.result = NULL,
                              path = NULL,
                              max.cols = 3,
                              plot.total = FALSE,
                              total.main = "System",
                              time.unit = "h",
                              layer.unit = "mg",
                              sink.unit = "ng",
                              mass.range = c(0, NULL),
                              time.range = c(0, NULL),
                              ...)
{
  # some input validation checks
  if (!.isInteger(max.cols) | max.cols < 1)
  {
    stop("max.col must be an integer and > 0.")
  }

  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  # read data
  if (is.null(r.result))
  {
    dat <- dskin.read.massData(parameter, path)
  }
  else
  {
    if (!.is_dskin_ds(r.result))
    {
      stop("r.result is not a DSkin data structure.")
    }

    dat <- r.result$masses
  }

  # Number of compartments
  n.data <- ncol(dat)
  if (n.data == 0)
  {
    stop("No logged data found. Cannot produce plots.")
  }

  if (plot.total)
  {
    n.data = n.data + 1
  }
  n.cols <- min(n.data, max.cols)
  n.rows <- ceiling(n.data/n.cols)

  mass.unit <- parameter$log$scaling
  old.par <- par(mfrow =c(n.rows, n.cols))

  total.mass <- rep(0, length(dat[[1]]$mass))
  total.time <- dat[[1]]$time

  sink.name <- parameter$compartments$sink$name
  x.lab = paste("time [", time.unit, "]", sep="")
  for (comp in names(dat))
  {
    values <- (dat[[comp]])
    time.vals <- dskin.convert.times(values$time, "min", time.unit)

    if (comp == sink.name)
    {
      mass.vals <- dskin.convert.masses(values$mass, mass.unit, sink.unit)
      fac = parameter$compartments$sink$Vd
      y.lab <- paste("c [", sink.unit, "/ml", "]", sep="")
    }
    else
    {
      fac = 1
      mass.vals <- dskin.convert.masses(values$mass, mass.unit, layer.unit)
      y.lab <- paste("mass [", layer.unit, "]", sep="")
    }

    x.range = range(time.vals)
    y.range = range(mass.vals)
    if (!is.null(time.range))
    {
      x.range[1] = ifelse(!is.na(time.range[1]), time.range[1], x.range[1])
      x.range[2] = ifelse(!is.na(time.range[2]), time.range[2], x.range[2])
    }
    if (!is.null(mass.range))
    {
      y.range[1] = ifelse(!is.na(mass.range[1]), mass.range[1], y.range[1])
      y.range[2] = ifelse(!is.na(mass.range[2]), mass.range[2], y.range[2])
    }

    plot(time.vals, mass.vals, main = comp, ..., xlab= x.lab, ylab=y.lab, xlim=x.range, ylim=y.range)

    total.mass <- total.mass + values$mass * fac
  }

  if (plot.total)
  {
    total.time <- dskin.convert.times(total.time, "min", time.unit)
    total.mass <- dskin.convert.masses(total.mass, mass.unit, layer.unit)

    x.range = range(total.time)
    y.range = range(total.mass)
    if (!is.null(time.range))
    {
      x.range[1] = ifelse(!is.na(time.range[1]), time.range[1], x.range[1])
      x.range[2] = ifelse(!is.na(time.range[2]), time.range[2], x.range[2])
    }
    if (!is.null(mass.range))
    {
      y.range[1] = ifelse(!is.na(mass.range[1]), mass.range[1], y.range[1])
      y.range[2] = ifelse(!is.na(mass.range[2]), mass.range[2], y.range[2])
    }

    y.lab <- paste("mass [", mass.unit, "]", sep="")
    plot(total.time, total.mass, main = total.main, ..., xlab= x.lab, ylab=y.lab,
         xlim=x.range, ylim=y.range)
  }

  par(old.par)
}


#' Overview Concentration-depth Profile Plots
#'
#' Produces an overview plot for concentration-depth-profiles (CDP).
#'
#' @param parameter The DSkin parameter pack.
#' @param r.result DSkin result data structure ('collect.data') by the R-version of DSkin. If NULL data is read from files defined by the Parameter Pack.
#' @param max.profiles The maximum number of profiles per plot (if not overridden by plotting.times).
#' @param plotting.times A vector of plotting times (time for every CDP). Overrides max.profiles if defined.
#' @param path The path to CDP files (NULL for wd).
#' @param max.cols The maximum number of columns in the plot.
#' @param layer.unit The mass unit for all layer compartments [ng,ug,mg].
#' @param conc.range The range of plottable concentration (NULL for automatic min/max).
#' @param space.range The space range for the plots (NULL for automatic min/max).
#' @param colors Vector of colors (or single color) for the color ramping procedure to colorize all curves in one subplot.
#' @param ... Additional parameter for plot directive.
#'
#' @return None
#'
#' @seealso \code{\link{dskin.parameter}}
#' @seealso \code{\link{dskin.simulate}}
#'
#' @importFrom graphics plot
#' @importFrom graphics par
#' @importFrom grDevices colorRampPalette
#' @export
dskin.plot.concs <- function(parameter,
                             r.result = NULL,
                             max.profiles = 5,
                             plotting.times = NULL,
                             path = NULL,
                             max.cols = 2,
                             layer.unit = "mg",
                             conc.range = c(0, NULL),
                             space.range = c(0, NULL),
                             colors=c("red", "yellow", "springgreen", "royalblue"),
                             ...)
{
  if (!.isInteger(max.cols) | max.cols < 1)
  {
    stop("max.col must be an integer and > 0.")
  }

  if (!.isInteger(max.profiles) | max.profiles < 1)
  {
    stop("max.profiles must be an integer and > 0.")
  }

  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  if (is.null(r.result))
  {
    dat <- dskin.read.concData(parameter, path)
  }
  else
  {
    if (!.is_dskin_ds(r.result))
    {
      stop("r.result is not a DSkin data structure.")
    }

    dat <- r.result$concs
  }

  n.data <- length(dat)
  if (n.data == 0)
  {
    stop("No logged data found. Cannot produce plots.")
  }

  n.cols <- min(n.data, max.cols)
  n.rows <- ceiling(n.data/n.cols)

  mass.unit <- parameter$log$scaling
  y.lab <- paste("c [", layer.unit, "/ml", "]", sep="")

  old.par <- par(mfrow =c(n.rows, n.cols), xpd=TRUE)
  # scaling data
  x.range = vector("list", n.data)
  y.range = vector("list", n.data)
  i <- 1
  for (comp in names(dat))
  {
    dat[[comp]]$data <- dskin.convert.masses(dat[[comp]]$data, mass.unit, layer.unit)
    x.r = range(dat[[comp]]$x)
    y.r = range(dat[[comp]]$data)
    if (!is.null(conc.range))
    {
      y.r[1] = ifelse(!is.na(conc.range[1]), conc.range[1], y.r[1])
      y.r[2] = ifelse(!is.na(conc.range[2]), conc.range[2], y.r[2])
    }
    if (!is.null(space.range))
    {
      x.r[1] = ifelse(!is.na(space.range[1]), space.range[1], x.r[1])
      x.r[2] = ifelse(!is.na(space.range[2]), space.range[2], x.r[2])
    }
    x.range[[i]] = x.r
    y.range[[i]] = y.r
    i <- i + 1
  }
  names(x.range) <- names(dat)
  names(y.range) <- names(dat)

  for (comp in names(dat))
  {
    col.fn <- colorRampPalette(colors)

    if (!is.null(plotting.times))
    {
      times <- plotting.times
    }
    else
    {
      times.length <- length(dat[[comp]]$time)
      times <- seq(0, times.length - 1, by=max(1, floor(times.length / max.profiles)))
    }
    x.vals <- dat[[comp]]$x
    y.vals <- dat[[comp]]$data[times + 1]

    ramp.colors <- col.fn(length(times))
    plot(x.vals,
         t(y.vals[1]),
         main = comp,
         ...,
         xlab= expression(paste("space [", mu, "m]", sep="")),
         ylab=y.lab,
         ylim=y.range[[comp]],
         xlim=x.range[[comp]],
         col=ramp.colors[1])

    if (length(times) > 1)
    {
      for (val in 2:length(times))
      {
        par(new = TRUE)
        plot(x.vals,
             t(y.vals[val]),
             ...,
             col=ramp.colors[val],
             ylim=y.range[[comp]],
             xlim=x.range[[comp]],
             xlab="",
             ylab="",
             axes=FALSE)
      }
    }
  }

  par(old.par)
}


#' Movie Creation for Mass-over-time profiles
#'
#' Produces an overview movie for all mass-over time profiles for all compartments.
#' This will produce temporary images in a subdirectory (relative to the current working directory).
#' Naming of the subdirectory follows the template 'TAG_mass_movie_tmp'. Here, TAG is pulled from the tag of the supplied
#' DSkin Parameter Pack. If the directory does exist it will not be deleted or cleared before the image creation started.
#'
#' @param parameter The DSkin parameter pack.
#' @param r.result DSkin result data structure ('collect.data') by the R-version of DSkin. If NULL data is read from files defined by the Parameter Pack.
#' @param remove.tmp.images Logical parameter to indicate whether produced temporary images should be deleted after processing.
#' @param path The path to the mass-over-time files (NULL for wd).
#' @param max.cols The maximum number of columns in the plot.
#' @param screen.play The screen play time of the produced movie in seconds.
#' @param frame.rate The appoximate frame rate of the produced movie.
#' @param resolution The resolution of the movie in c(width, height).
#' @param time.unit The unit of time [sec,min,h].
#' @param layer.unit The unit of the y-axis [ng,ug,mg].
#' @param sink.unit The mass unit for the sink compartment [ng,ug,mg].
#' @param time.label.pos The position of the time label in the plots (see [graphics::legend()]).
#' @param mass.range The range of plottable mass (NULL for automatic min/max).
#' @param time.range The time range for the plots (NULL for automatic min/max).
#' @param ... Additional parameter for plot directive.
#'
#' @return Filename of the produced movie.
#'
#' @seealso \code{\link{dskin.parameter}}
#' @seealso \code{\link{dskin.simulate}}
#'
#' @importFrom graphics plot
#' @importFrom graphics par
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom graphics legend
#' @importFrom stats window
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom doParallel stopImplicitCluster
#'
#' @export
dskin.movie.masses <- function(parameter,
                               r.result = NULL,
                               remove.tmp.images = TRUE,
                               path = NULL,
                               max.cols = 3,
                               screen.play = 10,
                               frame.rate = 25,
                               resolution = c(1280, 720),
                               time.unit = "h",
                               layer.unit = "mg",
                               sink.unit = "ng",
                               time.label.pos = "top",
                               mass.range = c(0, NULL),
                               time.range = c(0, NULL),
                               ...)
{
  if (!.isInteger(max.cols) | max.cols < 1)
  {
    stop("max.col must be an integer and > 0.")
  }

  if (!.dskin.check_in_path("ffmpeg"))
  {
    stop("Could not find ffmpeg in path.")
  }

  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  file.name <- paste(parameter$log$file_tag, "_masses.mp4", sep="")
  tmp.folder = paste(parameter$log$file_tag, "_mass_movie_tmp", sep="")

  # read data
  if (is.null(r.result))
  {
    dat <- dskin.read.massData(parameter, path)
  }
  else
  {
    if (!.is_dskin_ds(r.result))
    {
      stop("r.result is not a DSkin data structure.")
    }

    dat <- r.result$masses
  }

  n.data <- ncol(dat)
  if (n.data == 0)
  {
    stop("No logged data found. Cannot produce movie.")
  }

  base.dir <- getwd()
  if (!dir.exists(tmp.folder))
  {
    dir.create(tmp.folder, recursive = TRUE)
  }


  n.cols <- min(n.data, max.cols)
  n.rows <- ceiling(n.data/n.cols)
  mass.unit <- parameter$log$scaling

  max.imgs <- screen.play * frame.rate
  times <- dat[[1]]$time
  times.l <- length(times)
  if (times.l > max.imgs)
  {
    # thin out to maximize production speed
    # we must have the first time and last time on screen
    times <- window(times, deltat = floor(times.l/max.imgs))
    times.l <- length(times)
  }

  times.digits <- nchar(trunc(times.l))
  times.pf <- formatC(seq(1:times.l), width=times.digits, format="g", flag="0")
  setwd(tmp.folder)

  registerDoParallel(cores = detectCores())
  i <- 1 # to shut up devtools check
  foreach(i = 1:times.l) %dopar%
  {
    png(paste("mass_movie_", times.pf[i], ".png", sep=""),
        width=resolution[1], height = resolution[2])
    # write some time information
    old.par <- par(mfrow =c(n.rows, n.cols))
    sink.name <- parameter$compartments$sink$name
    for (comp in names(dat))
    {
      l.unit <- "mg"
      if (comp == sink.name)
      {
        l.unit = sink.unit
        y.lab <- paste("c [", sink.unit, "/ml", "]", sep="")
      }
      else
      {
        l.unit = layer.unit
        y.lab <- paste("mass [", layer.unit, "]", sep="")
      }

      values <- (dat[[comp]])

      x.range = dskin.convert.times(c(min(values$time), max(values$time)), "min", time.unit)
      y.range = dskin.convert.masses(c(min(values$mass), max(values$mass)), mass.unit, l.unit)
      if (!is.null(time.range))
      {
        x.range[1] = ifelse(!is.na(time.range[1]), time.range[1], x.range[1])
        x.range[2] = ifelse(!is.na(time.range[2]), time.range[2], x.range[2])
      }
      if (!is.null(mass.range))
      {
        y.range[1] = ifelse(!is.na(mass.range[1]), mass.range[1], y.range[1])
        y.range[2] = ifelse(!is.na(mass.range[2]), mass.range[2], y.range[2])
      }

      end.time <- times[i]
      end.time.idx <- which(values$time == end.time)

      time.data = dskin.convert.times(values$time[1:end.time.idx], "min", time.unit)
      mass.data = dskin.convert.masses(values$mass[1:end.time.idx], mass.unit, l.unit)

      plot(time.data,
           mass.data,
           main = comp,
           ...,
           xlab= paste("time [", time.unit, "]", sep=""),
           ylab= y.lab,
           xlim = x.range,
           ylim = y.range)

      t.text = paste(format(dskin.convert.times(times[i], "min", time.unit), nsmall = 2, digits=2),
                     time.unit)
      legend(time.label.pos, t.text, bty="n")
    }
    par(old.par)
    dev.off()
  }
  stopImplicitCluster()
  setwd(base.dir)

  file.pattern <- paste(tmp.folder, "/mass_movie_%0", times.digits, "d.png", sep="")
  sys.call <- "ffmpeg -framerate F_RATE -t SCREEN_PLAY -i PATTERN -c:v libx264 -profile:v high -crf 18 -pix_fmt yuv420p -y OUT_NAME"
  sys.call <- gsub("OUT_NAME", file.name, sys.call)
  sys.call <- gsub("PATTERN", file.pattern, sys.call)
  sys.call <- gsub("F_RATE", frame.rate, sys.call)
  sys.call <- gsub("SCREEN_PLAY", screen.play, sys.call)
  system(sys.call, ignore.stderr = TRUE)

  if (remove.tmp.images)
  {
    unlink(tmp.folder, recursive = TRUE)
  }

  return(file.name)
}


#' Movie Creation for Concentration-depth profiles
#'
#' Produces an overview movie for all concentration-depth profiles (CDPs) for all logged compartments.
#' This will produce temporary images in a subdirectory (relative to the current working directory).
#' Naming of the subdirectory follows the template 'TAG_conc_movie_tmp'. Here, TAG is pulled from the tag of the supplied
#' DSkin Parameter Pack. If the directory does exist it will not be deleted or cleared before the image creation started.
#'
#' @param parameter The DSkin parameter pack.
#' @param r.result DSkin result data structure ('collect.data') by the R-version of DSkin. If NULL data is read from files defined by the Parameter Pack.
#' @param remove.tmp.images Logical parameter to indicate whether produced temporary images should be deleted after processing.
#' @param path The path to CDP files (NULL for wd).
#' @param max.cols The maximum number of columns in the plot.
#' @param screen.play The screen play time of the produced movie in seconds.
#' @param frame.rate The appoximate frame rate of the produced movie.
#' @param resolution The resolution of the movie in c(width, height).
#' @param time.unit The unit of time [sec,min,h].
#' @param layer.unit The unit of the y-axis [ng,ug,mg].
#' @param time.label.pos The position of the time label in the plots (see [graphics::legend()]).
#' @param conc.range The range of plottable concentration (NULL for automatic min/max).
#' @param space.range The space range for the plots (NULL for automatic min/max).
#' @param auto.max.scale Logical, if y-axis maximum scaling should be adjusted for every frame
#' (every plot will display the current maximum concentation of the compartment).
#' @param ... Additional parameter for plot directive.
#'
#' @return Filename of the produced movie.
#'
#' @seealso \code{\link{dskin.parameter}}
#' @seealso \code{\link{dskin.simulate}}
#'
#' @importFrom graphics plot
#' @importFrom graphics par
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom graphics legend
#' @importFrom stats window
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom doParallel stopImplicitCluster
#'
#' @export
dskin.movie.concs <- function(parameter,
                              r.result = NULL,
                              remove.tmp.images = TRUE,
                              path = NULL,
                              max.cols = 3,
                              screen.play = 10,
                              frame.rate = 25,
                              resolution = c(1280, 720),
                              time.unit = "h",
                              layer.unit = "mg",
                              time.label.pos = "top",
                              conc.range = c(0, NULL),
                              space.range = c(0, NULL),
                              auto.max.scale = FALSE,
                              ...)
{
  if (!.isInteger(max.cols) | max.cols < 1)
  {
    stop("max.col must be an integer and > 0.")
  }

  if (!.dskin.check_in_path("ffmpeg"))
  {
    stop("Could not find ffmpeg in path.")
  }

  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  file.name <- paste(parameter$log$file_tag, "_concs.mp4", sep="")
  tmp.folder = paste(parameter$log$file_tag, "_conc_movie_tmp", sep="")

  # read data
  if (is.null(r.result))
  {
    dat <- dskin.read.concData(parameter, path)
  }
  else
  {
    if (!.is_dskin_ds(r.result))
    {
      stop("r.result is not a DSkin data structure.")
    }

    dat <- r.result$concs
  }

  n.data <- length(dat)
  if (n.data == 0)
  {
    stop("No logged data found. Cannot produce movie.")
  }

  base.dir <- getwd()
  if (!dir.exists(tmp.folder))
  {
    dir.create(tmp.folder, recursive = TRUE)
  }

  n.cols <- min(n.data, max.cols)
  n.rows <- ceiling(n.data/n.cols)
  mass.unit <- parameter$log$scaling

  max.imgs <- screen.play * frame.rate
  times <- dat[[1]]$time
  times.l <- length(times)
  if (times.l > max.imgs)
  {
    # thin out to maximize production speed
    # we must have the first time and last time on screen
    times <- window(times, deltat = floor(times.l/max.imgs))
    times.l <- length(times)
  }

  times.digits <- nchar(trunc(times.l))
  times.pf <- formatC(seq(1:times.l), width=times.digits, format="g", flag="0")
  setwd(tmp.folder)

  x.range = vector("list", n.data)
  y.range = vector("list", n.data)
  i <- 1
  for (comp in names(dat))
  {
    dat[[comp]]$data <- dskin.convert.masses(dat[[comp]]$data, mass.unit, layer.unit)
    x.r = range(dat[[comp]]$x)
    y.r = range(dat[[comp]]$data)
    if (!is.null(conc.range))
    {
      y.r[1] = ifelse(!is.na(conc.range[1]), conc.range[1], y.r[1])
      y.r[2] = ifelse(!is.na(conc.range[2]), conc.range[2], y.r[2])
    }
    if (!is.null(space.range))
    {
      x.r[1] = ifelse(!is.na(space.range[1]), space.range[1], x.r[1])
      x.r[2] = ifelse(!is.na(space.range[2]), space.range[2], x.r[2])
    }
    x.range[[i]] = x.r
    y.range[[i]] = y.r
    i <- i + 1
  }
  names(x.range) <- names(dat)
  names(y.range) <- names(dat)

  y.lab <- paste("c [", layer.unit, "/ml", "]", sep="")

  registerDoParallel(cores = detectCores())
  foreach(i=1:times.l) %dopar%
  {
    png(paste("conc_movie_", times.pf[i], ".png", sep=""),
        width=resolution[1], height = resolution[2])
    # write some time information
    old.par <- par(mfrow =c(n.rows, n.cols))
    for (comp in names(dat))
    {

      values <- dat[[comp]]

      end.time <- times[i]
      end.time.idx <- which(values$time == end.time)

      x.data <- values$x
      y.data <- values$data[end.time.idx]

      y.lim = y.range[[comp]]
      if (auto.max.scale)
      {
        y.lim[2] = max(y.data)
      }

      plot(x.data,
           t(y.data),
           main = comp,
           ...,
           xlab= expression(paste("space [", mu, "m]", sep="")),
           ylab= y.lab,
           xlim = x.range[[comp]],
           ylim = y.lim)

      t.text = paste(format(dskin.convert.times(values$time[end.time], "min", time.unit), nsmall = 2, digits=2),
                     time.unit)
      legend(time.label.pos, t.text, bty="n")
    }
    par(old.par)
    dev.off()
  }
  stopImplicitCluster()
  setwd(base.dir)

  file.pattern <- paste(tmp.folder, "/conc_movie_%0", times.digits, "d.png", sep="")
  sys.call <- "ffmpeg -framerate F_RATE -t SCREEN_PLAY -i PATTERN -c:v libx264 -profile:v high -crf 18 -pix_fmt yuv420p -y OUT_NAME"
  sys.call <- gsub("OUT_NAME", file.name, sys.call)
  sys.call <- gsub("PATTERN", file.pattern, sys.call)
  sys.call <- gsub("F_RATE", frame.rate, sys.call)
  sys.call <- gsub("SCREEN_PLAY", screen.play, sys.call)
  system(sys.call, ignore.stderr = TRUE)

  if (remove.tmp.images)
  {
    unlink(tmp.folder, recursive = TRUE)
  }

  return(file.name)
}


#' Converts a Movie to Gif
#'
#' Coverts a produced mp4 movie to gif file. The function needs ffmpeg in Path to work.
#' The output filename is constructed as 'TAG_masses.gif' for mass-over-time data or 'TAG_concs.gif' for concentration-depth profile data.
#' Here, TAG is the project tag defined in the DSkin Parameter Pack.
#'
#' @param parameter The DSkin Parameter Pack used to produce the source mp4 movie.
#' @param dskin.class The class of movie as a source [mass, conc].
#' @param width The desired width of the gif.
#' @param fps The desired FPS of the gif.
#'
#' @return Output filename of the produced gif.
#'
#' @seealso \code{\link{dskin.parameter}}
#'
#' @export
dskin.movie.to.gif <- function(parameter, dskin.class = "mass", width = 500, fps = 25)
{
  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin parameter pack.")
  }

  if (!.isInteger(width) | width < 100)
  {
    stop("width must be an integer and >= 100.")
  }

  if (!.isInteger(fps) | fps < 1)
  {
    stop("fps must be an integer and >= 1.")
  }

  if (!.dskin.check_in_path("ffmpeg"))
  {
    stop("Could not find ffmpeg in path.")
  }

  if (dskin.class == "mass")
  {
    source.filename <- paste(parameter$log$file_tag, "_masses.mp4", sep="")
    out.filename <- paste(parameter$log$file_tag, "_masses.gif", sep="")
    out.palette <- paste(parameter$log$file_tag, "_masses_palette.png", sep="")
  }
  else if (dskin.class == "conc")
  {
    source.filename <- paste(parameter$log$file_tag, "_concs.mp4", sep="")
    out.filename <- paste(parameter$log$file_tag, "_concs.gif", sep="")
    out.palette <- paste(parameter$log$file_tag, "_concs_palette.png", sep="")
  }
  else
  {
    stop("Unknown dskin.class. Allowed classes: mass, conc")
  }

  # first pass -> palette
  palette.call <- "ffmpeg -y -i IN_FILE -vf fps=FPS_VAL,scale=WIDTH:-1:flags=lanczos,palettegen PAL_OUT"
  palette.call <- gsub("IN_FILE", source.filename, palette.call)
  palette.call <- gsub("WIDTH", width, palette.call)
  palette.call <- gsub("PAL_OUT", out.palette, palette.call)
  palette.call <- gsub("FPS_VAL", fps, palette.call)
  system(palette.call, ignore.stderr = TRUE)

  # create gif
  gif.call <- "ffmpeg -y -i IN_FILE -i PAL_FILE -filter_complex \"fps=FPS_VAL,scale=WIDTH:-1:flags=lanczos[x];[x][1:v]paletteuse\" OUT_FILE"
  gif.call <- gsub("IN_FILE", source.filename, gif.call)
  gif.call <- gsub("PAL_FILE", out.palette, gif.call)
  gif.call <- gsub("WIDTH", width, gif.call)
  gif.call <- gsub("FPS_VAL", fps, gif.call)
  gif.call <- gsub("OUT_FILE", out.filename, gif.call)
  system(gif.call, ignore.stderr = TRUE)

  unlink(out.palette)

  return(out.filename)
}


#' Steady-state Line Producer
#'
#' Returns intercept and slope of a line to represent the steady-state of data produced by the DSkin steady-state finder.
#'
#' @param parameter The DSkin Parameter Pack.
#' @param ss.data DSkin steady-state data structure of class 'ss.fit.data'.
#' @param mass.unit Unit of mass of the y-axis in mass/ml. Possible values are [ng,ug,mg].
#' @param time.unit Unit of the x-axis [sec, min, h].
#'
#' @return Vector of interpect and slope or NA if no data could be produced.
#'
#' @seealso \code{\link{dskin.parameter}}
#' @seealso \code{\link{dskin.ss.analyze}}
#'
#' @export
dskin.ss.line <- function(parameter,
                          ss.data,
                          mass.unit = "mg",
                          time.unit = "h")
{
  if (!.is_dskin_ss_fit_data(ss.data))
  {
    stop("ss.data is not DSkin steady-state fit data.")
  }

  if (!.is_dskin_parameter(parameter))
  {
    stop("parameter is not a DSkin Parameter Pack.")
  }

  slope <- ss.data[row.names(ss.data) == "j.ss", ]
  if (is.na(slope))
  {
    return (c(NA, NA))
  }

  y.val <- 0
  x.val <- ss.data[row.names(ss.data) == "time.lag", ]

  ss.time.unit <- attr(ss.data, "time.unit")
  x.val <- dskin.convert.times(x.val, ss.time.unit, time.unit)

  ss.mass.unit <- attr(ss.data, "mass.unit")
  ss.space.unit <- attr(ss.data, "space.unit")
  # slope is in ss.mass.unit / (ss.time.unit * ss.space.unit**2)
  # mass in sink is always in unit / ml = unit / (cm^3)
  slope <- dskin.convert.masses(slope, ss.mass.unit, mass.unit)
  slope <- slope / dskin.convert.times(1, ss.time.unit, time.unit)

  # in cm^2
  A = parameter$compartments$vehicle$app_area
  # in ml
  V = parameter$compartments$sink$Vd
  slope <- slope / dskin.convert.spaces(1, ss.space.unit, "cm")**2
  slope <- slope * A
  slope <- slope / V
  intercept <- y.val - slope * x.val
  return (c(intercept, slope))
}



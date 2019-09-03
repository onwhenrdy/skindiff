library(devtools)
library(roxygen2)
library(Rcpp)

clean.up.src <- function()
{
  del.files <- list.files("./src", include.dirs = FALSE, full.names = TRUE)
  del.ignore <- c("RcppExports.cpp", "dskin_rbinding.cpp", "Makevars")
  del.files <- del.files[-grep(paste(del.ignore, collapse ="|"), del.files)]
  file.remove(del.files)
}

copy.src.files <- function()
{
  src.ignore <- c("main.cpp", "session.cpp", "session.h", "systemcmd.h", "systemcmd.cpp", "consoleprogressbar.h", "consoleprogressbar.cpp")
  src.include <- c(".cpp", ".h")
  src.files <- list.files("../../", include.dirs = FALSE, full.names = TRUE)
  src.files <- src.files[grep(paste(src.include, collapse ="|"), src.files)]
  src.files <- src.files[-grep(paste(src.ignore, collapse ="|"), src.files)]
  file.copy(src.files, "./src", overwrite =TRUE)
}


####
base.dir <- getwd()
setwd("./dskin")
copy.src.files()

if (.Platform$OS.type == "windows")
{
  find_rtools()
}

pkgbuild::compile_dll()
check.res <- check(quiet = FALSE)
if (!is.null(check.res$errors) && length(check.res$errors) > 0)
{
  stop(check.res$errors)
}

if (!is.null(check.res$warnings) && length(check.res$warnings) > 0)
{
  # we only allow "non-portable flags" warning
  if (length(check.res$warnings) == 1 && length(grep("Non-portable flags", check.res$warnings)) < 1)
  {
    stop(check.res$warnings)
  }
}

if (!is.null(check.res$notes) && length(check.res$notes) > 0)
{
  stop(check.res$notes)
}

pkg.file <- build(vignettes = FALSE, quiet = TRUE)
cat(paste("Build Package: ", pkg.file, "\n"))

clean.up.src()
setwd(base.dir)

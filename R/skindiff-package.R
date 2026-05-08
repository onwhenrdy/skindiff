#' skindiff: numerical simulation of drug diffusion through skin
#'
#' Builds a 1-D Crank-Nicolson finite-difference model of a stack of skin
#' compartments (vehicle / stratum corneum / deeper skin layers / sink) and
#' returns time-resolved mass and concentration-depth profiles.
#'
#' @useDynLib skindiff, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

#' @keywords internal
.dskin.sys_params <- function()
{
  tmp = list("disc_scheme" = "EQUIDIST",
             "mb_method" = "DSkin_1_4",
             "resolution" = 1,
             "max_module" = 50.0,
             "mb_eta" = 0.6,
             "sim_time" = 600)
  attr(tmp, "name") <- "sys"

  return(tmp)
}


#' @keywords internal
.dskin.pk_params <- function()
{
  tmp = list("enabled" = TRUE,
             "t_half" = 1.0)
  attr(tmp, "name") <- "PK"

  return(tmp)
}


#' @keywords internal
.dskin.log_params <- function(file_tag)
{
  f_tag = "test"
  if (!missing(file_tag))
  {
    if (!.isSingleString(file_tag))
    {
      stop("file_tag is not a single string.")
    }

    f_tag = file_tag
  }

  tmp = list("file_tag" = f_tag,
             "mass_file_postfix" = "mass",
             "mass_file_gzip" = FALSE,
             "cdp_file_postfix" = "cdp",
             "cdp_file_gzip" = TRUE,
             "mass_log_interval" = 1,
             "cdp_log_interval" = 1,
             "scaling" = "mg",
             "show_progress" = TRUE)
  attr(tmp, "name") <- "log"

  return(tmp)
}


#' @keywords internal
.dskin.vehicle_params <- function(comp_name)
{
  c_name = "Donor"
  if (!missing(comp_name))
  {
    if (!.isSingleString(comp_name))
    {
      stop("comp_name is not a single string.")
    }

    c_name = comp_name
  }

  tmp = list("name" = c_name,
             "finite_dose" = TRUE,
             "c_init" = 1.0,
             "app_area" = 1.0,
             "h" = 30,
             "D" = 1.0,
             "replace_after" = 0,
             "remove_after" = 0,
             "log" = TRUE,
             "log_cdp" = TRUE)
  attr(tmp, "name") <- "vehicle"

  return(tmp)
}


#' @keywords internal
.dskin.sink_params <- function(comp_name)
{
  c_name = "Sink"
  if (!missing(comp_name))
  {
    if (!.isSingleString(comp_name))
    {
      stop("comp_name is not a single string.")
    }

    c_name = comp_name
  }

  tmp = list("name" = c_name,
             "log" = TRUE,
             "c_init" = 0.0,
             "Vd" = 1.0)
  attr(tmp, "name") <- "sink"

  return(tmp)
}


# compile parameter packs to a consistent list form
#' @keywords internal
.dskin.compile_params <- function(...)
{
  tmp.list <- list(...)
  tmp.list[sapply(tmp.list, is.null)] <- NULL
  names(tmp.list) <- sapply(tmp.list, function(x) {attr(x, "name")})

  return(tmp.list)
}


#' @keywords internal
.dskin.layer_params <- function(n_layers, layer_names)
{
  n_lay = 0
  if (!missing(n_layers))
  {
    if(!.isInteger(n_layers) | n_layers < 0)
    {
      stop("n_layers is not an integer >= 0.")
    }

    n_lay = n_layers;
  }

  l_names <- c()
  has_layer_names = missing(layer_names) | is.null(layer_names)
  if (n_lay > 0 & has_layer_names)
  {
    l_names <- sprintf("Layer%d",seq(1:n_lay))
  }

  if (n_lay > 0 & !has_layer_names)
  {
    if (length(layer_names) != n_lay)
    {
      stop("Number of layers != length of layer_names")
    }

    if (!all(sapply(layer_names, .isSingleString)))
    {
      stop("layer_names must be a list of strings.")
    }

    if (anyDuplicated(layer_names))
    {
      stop("layer_names must be unique. Found duplicates.")
    }

    l_names = layer_names
  }

  if (n_lay == 0 & !has_layer_names)
  {
    if (length(layer_names) != 0)
    {
      stop("0 layers defined but non-empty layer_names found.")
    }
  }

  if (n_lay != 0)
  {
    l_vec = c()
    for (i in seq(1, n_lay))
    {
      layer = list("name" = l_names[i],
                   "log" = TRUE,
                   "log_cdp" = TRUE,
                   "c_init" = 0.0,
                   "cross_section" = 1.0,
                   "h" = 10,
                   "D" = 1.0,
                   "K" = 1.0)

      l_vec = c(l_vec, layer)
    }

    names.list <- c("name", "log", "log_cdp", "c_init", "cross_section", "h", "D", "K")
    df <- data.frame(matrix(unlist(l_vec), nrow=n_lay, byrow=T), stringsAsFactors=FALSE)
    colnames(df) <- names.list
    df[, c(2,3)] <- sapply(df[, c(2,3)], as.logical)
    df[, c(4,5,7,8)] <- sapply(df[, c(4,5,7,8)], as.double)
    df[, c(6)] <- sapply(df[, c(6)], as.integer)

    attr(df, "name") <- "layers"
    return (df)
  }
}


#' @keywords internal
.dskin.compartment_params <- function(n_layers, layer_names, vehicle_name, sink_name)
{
  v_name <- "Donor"
  if (!missing(vehicle_name))
  {
    v_name <- vehicle_name
  }

  s_name <- "Sink"
  if (!missing(sink_name))
  {
    s_name <- sink_name
  }

  vehicle = .dskin.vehicle_params(v_name)
  sink = .dskin.sink_params(s_name)
  layers = .dskin.layer_params(n_layers, layer_names)

  comps <- .dskin.compile_params(vehicle, sink, layers)
  attr(comps, "name") <- "compartments"

  return (comps)
}

#' Calculates the temperature sensitivity of a model (and parameters)
#'
#' @param data a nested list of data formatted using one of the format_*()
#' functions
#' @param par a vector of parameter values, this is functions specific
#' @param model the model name to be used in optimizing the model
#' @param ... extra arguments to pass to the function
#' @return returns a temperature sensitivity as days per degree celsius
#' @keywords phenology, model, post-processing
#' @export
#' @examples
#'
#' \dontrun{
#'
#' # model calibration with default parameters
#' # phenocam_DB, "TT" to estimate a parameter set
#' pars <- pr_fit()
#'
#' # estimate temperatue sensitivity for a particular
#' # combination of data, parameters and model
#' sensitivity <- pr_calc_temp_sens(par = pars$par,
#'  data = phenocam_DB, model = "TT")
#'
#' }

pr_calc_temp_sens <- function(
  par,
  data,
  model,
  ...
  ){

  # check parameters
  if(missing(par)){
    stop("please provide a parameter set")
  }

  # check data
  if(missing(data)){
    stop("please provide a valid data file or path")
  }

  # check model
  if(missing(model)){
    stop("please provide a model name")
  }

  # when provided with a file name,
  # load the rds file (name data)
  if(is.character(data)){
    data <- readRDS(data)
  }

  # convert to a flat format (to be sure) and duplicate
  data_t1 <- data_t0 <- flat_format(data = data)

  # increase all temperatur related paramters with 1 degree C
  data_t1$Ti <- data_t1$Ti + 1
  data_t1$Tmini <- data_t1$Tmini + 1
  data_t1$Tmaxi <- data_t1$Tmaxi + 1

  # return model output with normal data (t0)
  t0 <- do.call(model, list(data = data_t0,
                        par = par,
                        ... ))

  # return model output with warmer data (t1)
  t1 <- do.call(model, list(data = data_t1,
                            par = par,
                            ... ))

  # set out of range values to NA
  t0[t0 == 9999] <- NA
  t1[t1 == 9999] <- NA

  # calculate difference in doy when increasing
  # temperatures with 1 degree C
  doy_diff <- t1 - t0

  # return the raw data and summary stats for
  # convenience
  return(list(mean = mean(doy_diff, na.rm = TRUE),
              sd = stats::sd(doy_diff, na.rm = TRUE),
              doy_diff = doy_diff))
}

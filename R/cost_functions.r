#' Root Mean Squared Error cost function for model optimization.
#'
#' @param par a vector of parameter values, this is functions specific
#' @param data nested data structure with validation data as returned
#' by format_phenocam() or format_pep725(), or your own dataset adhering
#' to the same data structure.
#' @param model the model name to be used in optimizing the model
#' @param ... extra arguments to pass to the function
#' @return the RMSE comparing observed and estimated values
#' @keywords phenology, model, optimization, cost function
#' @export
#' @examples
#'
#' # The cost function returns the rmse between the
#' # true values and those generated by the model given a
#' # parameterset par.
#' \dontrun{
#' cost_value = rmse(par, data, model="TTs")
#' }

rmse <- function(
  par,
  data,
  model,
  ...
  ) {

  # inset validity checks
  val <- data$transition_dates
  out <- do.call(model, list(data = data, par = par, ...))

  if (any(is.na(out))) {
    return(9999)
  } else {
    # return the RMSE between the validation data and
    # the output of the model
    return(sqrt(mean((val - out) ^ 2, na.rm = T)))
  }
}

#' Coefficient of Variation of the Mean Absolute Error
#'
#' @param par a vector of parameter values, this is functions specific
#' @param data nested data structure with validation data as returned
#' by format_phenocam() or format_pep725(), or your own dataset adhering
#' to the same data structure.
#' @param model the model name to be used in optimizing the model
#' @param ... extra arguments to pass to the function
#' @return the CVMAE comparing observed and estimated values
#' @keywords phenology, model, optimization, cost function
#' @export
#' @examples
#'
#' # The cost function returns the CVMAE between the
#' # true values and those generated by the model given a
#' # parameterset par.
#' \dontrun{
#' cost_value = cvmae(par, data, model="TTs")
#' }

cvmae <- function(
  par,
  data,
  model,
  ...
  ) {

  # combine in data frame
  df <- data.frame(
    val = data$transition_dates,
    out = do.call(model, list(data = data, par = par, ...)),
    site = data$site
  )

  if (any(is.na(df$out))) {
    return(9999)
  } else {
    return(
      mean(by(df, INDICES = df$site, function(x){
        mean(abs(x$val - x$out))/mean(x$val)
      }, simplify = TRUE))
    )
  }
}

#' Log likelihood cost function for model optimization
#'
#' The function is aimed to be maximized, to use it with optimizers which
#' minimize cost functions wrap the function as such:
#' `cost = function(...)\{abs(likelihood(...))\}`
#'
#' @param par a vector of parameter values, this is functions specific
#' @param data nested data structure with validation data as returned
#' by format_phenocam() or format_pep725(), or your own dataset adhering
#' to the same data structure.
#' @param model the model name to be used in optimizing the model
#' @param ... extra arguments to pass to the function
#' @return the RMSE comparing observed and estimated values
#' @keywords phenology, model, optimization, cost function
#' @export
#' @examples
#'
#' # The cost function returns the rmse between the
#' # true values and those generated by the model given a
#' # parameterset par.
#' \dontrun{
#' cost_value = likelihood(par, data, model="TTs")
#' }

likelihood <- function(
  par,
  data,
  model,
  ...
  ) {

  # model parameters
  model_par <- par[1:(length(par)) - 1]

  # split out sd range parameter
  sd_range <- par[length(par)]

  # run model
  observed <- data$transition_dates
  predicted <- do.call(
    model,
    list(
      data = data,
      par = model_par
      )
    )

  # get residuals
  residuals <- predicted - observed

  # singlelikelihood
  singlelikelihoods <- stats::dnorm(
      residuals,
      sd = sd_range,
      log = TRUE
      )

  return(sum(singlelikelihoods, na.rm = TRUE))
}

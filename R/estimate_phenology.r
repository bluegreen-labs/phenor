#' Calculates the estimated phenophases for a given set of parameters
#' and a specified model (be sure to match parameter and requirements
#' with the model. Wrapper around a do.call() call for convenience.
#'
#' @param data: a nested list of data with on location:
#' 1. the date (doy or long format)
#' 2. the temperature data
#' 3. the photoperiod data (NA when not needed)
#' 4. a vector or matrix with necessary constants (NA when not needed)
#'    - long term mean temperature
#'    - latitude
#'    - etc...
#' 5. validation data (optional when just running a model not optimizing)
#' @param par: a vector of parameter values, this is functions specific
#' @param model: the model name to be used in optimizing the model
#' @keywords phenology, model, post-processing
#' @export
#' @examples
#' estimate <- estimate.phenology(par,data,model)
#'
#' # estimate will be an estimated timing of a phenophase

# calls a particular model and executes it using a given
# set of parameters and data. Make sure to check the data
# and parameter constraints. Data matrices must match.

estimate_phenology = function(par, data, model) {

  # convert to a flat format for speed
  data = flat_format(data = data)

  # return model output
  do.call(model, list(data = data, par = par))
}

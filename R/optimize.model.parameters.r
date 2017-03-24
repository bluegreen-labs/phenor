#' Calculates the estimated phenophases for a given set of parameters
#' and a specified model (be sure to match parameter and requirements
#' with the model.
#'
#' @param par: a vector of starting parameter values (function specific)
#' @param data: a nested list of data with on location:
#' 1. the date (doy or long format)
#' 2. the temperature data
#' 3. the photoperiod data (NA when not needed)
#' 4. a vector or matrix with necessary constants (NA when not needed)
#'    - long term mean temperature
#'    - latitude
#'    - etc...
#' 5. validation data (REQUIRED)
#' @param cost: the cost function to use in the optimization, it should return
#' a RMSE or other value which needs to be minimized
#' @param model: the model name to be used in optimizing the model
#' @param method: optimization method to use
#'    - GenSA :  Generalized Simulated Annealing algorithm (default)
#'    - genoud : GENetic Optimization Using Derivatives
#'    - those available in optim()
#'
#' @param lower: lower limit of parameter values (function specific)
#' @param upper: upper limit of parameter values (function specific)
#' @param maxit: maximum number of iterations to run (if NULL until convergence)
#' @keywords phenology, model, optimization, simulated annealing, genoud, optim
#' @export
#' @examples
#' estimate <- estimate.phenology(par,data,model)
#'
#' # estimate will return the best estimated parameter set given the
#' # validation data

optimize.model.parameters = function(par = par,
                                     data = data,
                                     cost = cost.function,
                                     model = "SEQ1.3",
                                     method = "GenSA",
                                     lower = -Inf,
                                     upper = Inf,
                                     maxit = 100){

  # load required libraries
  require(likelihood, quietly = TRUE)
  require(GenSA, quietly = TRUE)

  # solid results
  if ( method == "GenSA" ){
    # optimize model parameters using the
    # Generalized Simulated Annealing algorithm
    optim.par = GenSA(par = par,
                      fn = cost,
                      data = data,
                      lower = lower,
                      upper = upper,
                      control = list(maxit=maxit))
  }

  # does not work yet debug
  if (method == "genoud"){
    # optimize model parameters using the
    # GENetic Optimization Using Derivatives
    optim.par = genoud(fn = cost,
                       nvars = length(par),
                       data = data,
                       par = par,
                       max.generations = maxit)
  }

  # return the optimization data (parameters)
  return(optim.par)
}

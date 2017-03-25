#' Calculates the estimated phenophases for a given set of parameters
#' and a specified model (be sure to match parameter and requirements
#' with the model.
#'
#' @param par: a vector of starting parameter values (function specific)
#' @param data: a nested list of data with on location:
#' 1. the date (doy or long format)
#' 2. the long term daily mean temperature data
#' 3. the geographic location (lat, lon)
#' 4. a vector or matrix with necessary constants (NA when not needed)
#'    - long term mean temperature
#'    - latitude
#'    - etc...
#' 5. phenophase validation data (REQUIRED)
#' @param cost: the cost function to use in the optimization, it should return
#' a RMSE or other value which needs to be minimized
#' @param model: the model name to be used in optimizing the model
#' @param method: optimization method to use
#'    - GenSA :  Generalized Simulated Annealing algorithm (default)
#'    - genoud : GENetic Optimization Using Derivatives
#'    - those available in optim()
#' @param lower: lower limit of parameter values (function specific)
#' @param upper: upper limit of parameter values (function specific)
#' @param maxit: maximum number of iterations to run (if NULL until convergence)
#' @keywords phenology, model, optimization, simulated annealing, genoud, optim
#' @export
#' @examples
#'
#' estimate <- estimate.phenology(par,data,model)
#'
#' # estimate will return the best estimated parameter set given the
#' # validation data

optimize.parameters = function(par = NULL,
                               data = data,
                               cost = cost.function,
                               model = "SEQ1.3",
                               method = "GenSA",
                               lower = NULL,
                               upper = NULL,
                               maxit = 1000) {

  # load required libraries
  require(likelihood, quietly = TRUE)
  require(GenSA, quietly = TRUE)

  # check if starting parameters are present
  if (is.null(lower) | is.null(upper)){
    stop('Please provide upper and lower boundaries to the parameter space.\n
         Not defining your parameters space might yield good fits,\n
         for a wrong reason.')
  }

  if ( method == "GenSA" ){
    # one can opt to automatically generate starting values
    # in GenSA, this might yield better results. Set the
    # par parameter to NULL to do so.

    optim.par = GenSA(par = par,
                      fn = cost,
                      data = data,
                      lower = lower,
                      upper = upper,
                      control = list(maxit=maxit))
  }

  if (method == "genoud"){
    # optimize model parameters using the
    # GENetic Optimization Using Derivatives
    # needs more tweaking to work out of the box
    # on most models

    # stop if no starting parameters are provided
    if (is.null(par)){
      stop('The genoud algorithm needs defined strating parameters!')
    }

    optim.par = genoud(fn = cost,
                       nvars = length(par),
                       data = data,
                       par = par,
                       max.generations = maxit,
                       data.type.int = FALSE)
  }

  # other optimizers can be added here !

  # return the optimization data (parameters)
  # check formatting for post-processing
  return(optim.par)
}

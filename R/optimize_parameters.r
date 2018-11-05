#' Calculates the estimated phenophases for a given set of parameters
#' and a specified model (be sure to match parameter and requirements
#' with the model.
#'
#' @param par a vector of starting parameter values (function specific)
#' @param data nested data structure with validation data as returned
#' by format_phenocam() or format_pep725(), or your own dataset adhering
#' to the same data structure.
#' @param cost the cost function to use in the optimization, it should return
#' a RMSE or other value which needs to be minimized
#' @param model the model name to be used in optimizing the model
#' @param method optimization method to use (default = GenSA)
#'    - GenSA :  Generalized Simulated Annealing algorithm
#'    - genoud : GENetic Optimization Using Derivatives
#'    - BayesianTools: various bayesian based optimization tools
#' @param lower lower limit of parameter values (function specific)
#' @param upper upper limit of parameter values (function specific)
#' @param maxit maximum number of generations to run (genoud)
#' @param control additional optimization control parameters (default = NULL)
#' @param ... extra arguments to pass to the function, mostly BayesianTools
#' @keywords phenology, model, optimization, simulated annealing, genoud, optim
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- estimate.phenology(par,data,model)
#'
#' # estimate will return the best estimated parameter set given the
#' # validation data
#' }

optimize_parameters = function(par = NULL,
                               data = data,
                               cost = rmse,
                               model = "TT",
                               method = "GenSA",
                               lower = NULL,
                               upper = NULL,
                               maxit = NULL,
                               control = NULL,
                               ... ) {

  # check if starting parameters are present
  if (is.null(lower) | is.null(upper)){
    stop('Please provide upper and lower boundaries to the parameter space.\n
         Not defining your parameters space might yield good fits,\n
         for a wrong reason.')
  }

  # check if starting parameters are balanced
  if (length(lower) != length(upper)){
    stop('Parameter boundaries should be balanced')
  }

  # convert to a flat format for speed
  # (increases speed > 20x)
  data = flat_format(data)

  if ( tolower(method) == "gensa" ){
    # one can opt to automatically generate starting values
    # in GenSA, this might yield better results. Set the
    # par parameter to NULL to do so.
    optim_par = GenSA::GenSA(
      par = par,
      data = data,
      fn = cost,
      lower = lower,
      upper = upper,
      model = model,
      control = control,
      ...
    )
  }

  if (tolower(method) == "genoud"){
    # optimize model parameters using the
    # GENetic Optimization Using Derivatives
    # needs more tweaking to work out of the box
    # on most models

    # stop if no starting parameters are provided
    if (is.null(par)){
      stop('The genoud algorithm needs defined strating parameters!')
    }
    optim_par = rgenoud::genoud(
      fn = cost,
      nvars = length(par),
      max.generations = maxit,
      Domains = cbind(lower,upper),
      boundary.enforcement = 2,
      data.type.int = FALSE,
      data = data,
      model = model,
      ...
    )
  }

  # BayesianTools
  if (tolower(method) == "bayesiantools"){

    # Set the sd metric fixed to 1/5 of the range defined
    # by upper and lower limits. This ensures proper sampling
    # across widely different ranges.
    sd_range = abs(upper - lower)/5

    # setup the bayes run, no message forwarding is provided
    # so wrap the function in a do.call
    setup = BayesianTools::createBayesianSetup(
      likelihood = function(random_par){
        do.call("likelihood",
                list(par = random_par,
                     data = data,
                     model = model,
                     sd_range = sd_range))},
        lower = lower,
        upper = upper,
        ...
      )

      # calculate the runs
      out = BayesianTools::runMCMC(bayesianSetup = setup,
                                   sampler = control$sampler,
                                   settings = control$settings)

      # correct formatting in line with other outputs
      optim_par = list("par" = BayesianTools::MAP(out)$parametersMAP,
                       "bt_output" = out)
    }

  # return the optimization data (parameters)
  # check formatting for post-processing
  return(optim_par)
}

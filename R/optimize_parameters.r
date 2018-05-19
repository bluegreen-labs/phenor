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
#' @param method optimization method to use
#'    - GenSA :  Generalized Simulated Annealing algorithm (default)
#'    - genoud : GENetic Optimization Using Derivatives
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

  if ( method == "GenSA" ){
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

  if (method == "genoud"){
    # optimize model parameters using the
    # GENetic Optimization Using Derivatives
    # needs more tweaking to work out of the box
    # on most models

    # stop if no starting parameters are provided
    if (is.null(par)){
      stop('The genoud algorithm needs defined strating parameters!')
    }
    optim_par = rgenoud::genoud(fn = cost,
                       nvars = length(par),
                       max.generations = maxit,
                       Domains = cbind(lower,upper),
                       boundary.enforcement = 2,
                       data.type.int = FALSE,
                       data = data,
                       model = model,
                       ...)
  }

  # BayesianTools
  if (tolower(method) == "bayesiantools"){

    # dangerous shit
    assign("tmp_data", data, envir = .GlobalEnv)
    assign("tmp_model", model, envir = .GlobalEnv)

    # setup the bayes run
    setup = BayesianTools::createBayesianSetup(likelihood = likelihood,
                                 lower = c(lower, 0.01),
                                 upper = c(upper, 30))

    # run the MCMC routine, pass the sampler
    # via ...
    out = BayesianTools::runMCMC(bayesianSetup = setup,
                    settings = control)

    # remove copy of data and model name
    try(rm("tmp_data","tmp_model"))

    # correct formatting in line with other outputs
    optim_par = list("par" = BayesianTools::MAP(out)$parametersMAP[1:length(lower)])
  }

  # return the optimization data (parameters)
  # check formatting for post-processing
  return(optim_par)
}

#' Model comparison routine to faciliate model development
#' and quick comparisons of the skill of various models.
#'
#' @param random_seeds a vector with random seeds for cross validation
#' @param models list of models to compare
#' @param data which standard or custom dataset to use
#' @param method optimization method to use (default = GenSA)
#'    - GenSA :  Generalized Simulated Annealing algorithm
#'    - genoud : GENetic Optimization Using Derivatives
#'    - BayesianTools: various bayesian based optimization tools
#' @param control additional optimization control parameters
#' (default = list(max.call = 5000, temperature = 10000))
#' @param par_ranges location of the parameter ranges of the models
#' @param ncores number of cores to use to calculate model comparisons,
#' system specific and defaults to 1 threat (default = 1)
#' @keywords phenology, model, validation, comparison
#' @export
#' @examples
#'
#' # estimate will return the best estimated parameter set given the
#' # validation data
#' \dontrun{
#' my_comparison <- pr_fit_compare(random_seeds = c(38,1),
#'  models = c("TT","PTT"),
#'  dataset = "phenocam_DB",
#'  par_ranges = "parameter_ranges.csv")
#' }

pr_fit_compare <- function(
  random_seeds = c(1,12,40),
  models = c("LIN","TT","TTs","PTT","PTTs",
             "M1","M1s","AT","SQ","SQb","SM1",
             "SM1b","PA","PAb","PM1",
             "PM1b","UN","UM1","SGSI","AGSI"),
  data = phenor::phenocam_DB,
  method = "GenSA",
  control = list(max.call = 5000,
                 temperature = 10000),
  par_ranges = system.file(
    "extdata",
    "parameter_ranges.csv",
    package = "phenor",
    mustWork = TRUE
  ),
  ncores = 1
){

  # convert to a flat format for speed
  data = flat_format(data)

  # load parameter ranges
  par_ranges = utils::read.table(par_ranges,
                          header = TRUE,
                          sep = ",")

  # subset the parameter range
  if (!any(par_ranges$model %in% models)){
    stop("parameters are not specified in the standard ")
  }

  # get dimensions
  nr_models = length(models)
  nr_seeds = length(random_seeds)
  nr_rows = nr_models * nr_seeds
  nr_obs = length(data$transition_dates)

  # initiate list
  model_estimates = list()

  # Both for loops can be parallelized,
  # to speed up comparisons.

  if (tolower(method) != "bayesiantools"){

    # implement a progress bar for graphical feedback
    # this to gauge speed limitations
    cat("This might take a while ... \n")
    pb = utils::txtProgressBar(min = 0, max = nr_models, style = 3)
    k = 1

    # iterate all instances
    for (i in 1:nr_models) {

      # select ranges
      d = par_ranges[par_ranges$model == models[i],]
      d = d[,!is.na(d[1,])]
      d = d[,3:ncol(d)]
      d = as.matrix(d)

      # start cluster, default is SOCK cluster
      cl <- snow::makeCluster(ncores)

      # optimize models
      optimized_data = snow::parLapply(cl,
                                       random_seeds,
                                       function(random_seed){

        # load library in individual workers
        # similar the foreach .package argument
        library(phenor)

        # set random seed for a given run
        set.seed(random_seed)

        # progress bar for the models
        utils::setTxtProgressBar(pb, k);
        k = k + 1

        # optimize the model parameters using
        # GenSA algorithm
        par = pr_fit_parameters(
          par = NULL,
          data = data,
          model = models[i],
          method = method,
          lower = as.numeric(d[1,]),
          upper = as.numeric(d[2,]),
          control = control
        )

        # add model output using the estiamted parameters
        predicted_values = pr_predict(
          data = data,
          model = models[i],
          par = par$par
        )

        # return data
        return(list("parameters" = par$par,
                    "predicted_values" = predicted_values))
      })

      # stop cluster
      snow::stopCluster(cl)

      # stuff everything in a list
      tmp = list("parameters" = do.call("rbind",
                                        lapply(optimized_data,
                                               function(x)x$parameters)),
                 "predicted_values" = do.call("rbind",
                                              lapply(optimized_data,
                                                     function(x)x$predicted_values)))

      # append to the list
      model_estimates[i] = list(tmp)
    }

    # close the progress bar element
    close(pb)

  } else {

    # warning on random seeds
    if(length(random_seeds) > 1){
      message("Only the first random seed will be used,
              please use the 'nchains' parameter
              in BT to specify the number random initial conditions!")
    }

    # parallel process the models optimizations
    # set ncores to number of models if specified
    # ncores exceeds the number of models
    if(ncores > nr_models) ncores <- nr_models

    # implement a progress bar for graphical feedback
    # this to gauge speed limitations
    cat("This might take a while ... \n")

    # start cluster, default is SOCK cluster
    cl <- snow::makeCluster(ncores)

    # optimize models
    model_estimates = snow::parLapply(cl, models, function(model){

      # select ranges
      d = par_ranges[par_ranges$model == model,]
      d = d[,!is.na(d[1,])]
      d = d[,3:ncol(d)]
      d = as.matrix(d)

      # optimize the model parameters using
      # GenSA algorithm
      par = phenor::optimize_parameters(
        par = NULL,
        data = data,
        model = model,
        method = method,
        lower = as.numeric(d[1,]),
        upper = as.numeric(d[2,]),
        control = control
      )

      # add model output using the estiamted parameters
      predicted_values = estimate_phenology(
        data = data,
        model = model,
        par = par$par
      )

      # return data
      return(list("parameters" = par$par,
                  "predicted_values" = predicted_values,
                  "parameter_uncertainty" = par$par,
                  "bt_output" = par$bt_output))
    })

    # stop cluster
    snow::stopCluster(cl)

  }

  # collect garbage, especially unclosed connections
  # or memory stacks
  gc()

  # rename model output for clarity
  names(model_estimates) = models

  # concat comparison data
  comparison = list("modelled" = model_estimates,
                    "measured" = data$transition_dates)

  # return data
  return(comparison)
}

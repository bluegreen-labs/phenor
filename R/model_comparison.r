#' Model comparison routine to faciliate model development
#' and quick comparisons of the skill of various models.
#'
#' @param random_seeds a vector with random seeds for cross validation
#' @param models list of models to compare
#' @param dataset which standard or custom dataset to use
#' @param method which optimization method to use, GenSA or rgenoud
#' (default = GenSA)
#' @param control additional optimization control parameters
#' (default = list(max.call = 5000, temperature = 10000))
#' @param par_ranges location of the parameter ranges of the models
#' @keywords phenology, model, validation, comparison
#' @export
#' @examples
#'
#' # estimate will return the best estimated parameter set given the
#' # validation data
#' \dontrun{
#' my_comparison <- model_comparison(random_seeds = c(38,1),
#'  models = c("TT","PTT"),
#'  dataset = "phenocam_DB",
#'  par_ranges = "parameter_ranges.csv")
#' }

model_comparison = function(random_seeds = c(1,12,40),
                            models = c("LIN","TT","TTs","PTT","PTTs",
                                       "M1","M1s","AT","SQ","SQb","SM1",
                                       "SM1b","PA","PAb","PM1",
                                       "PM1b","UN","UM1","SGSI","AGSI"),
                            dataset = phenocam_DB,
                            method = "GenSA",
                            control = list(max.call = 5000,
                                           temperature = 10000),
                            par_ranges = sprintf("%s/extdata/parameter_ranges.csv",
                                                 path.package("phenor"))){

  # convert to a flat format for speed
  data = flat_format(dataset)

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

  # implement a progress bar for graphical feedback
  # this to gauge speed limitations
  cat("This might take a while ... \n")
  pb = utils::txtProgressBar(min = 0, max = nr_models*nr_seeds, style = 3)
  k = 1

  # iterate all instances
  for (i in 1:nr_models) {

    # select ranges
    d = par_ranges[par_ranges$model == models[i],]
    d = d[,!is.na(d[1,])]
    d = d[,3:ncol(d)]
    d = as.matrix(d)

    # create temp matrices
    tmp_parameters = matrix(NA,nr_seeds,ncol(d))
    predicted_values = matrix(NA,nr_seeds,nr_obs)

    for (j in 1:nr_seeds) {

      # progress bar for the models
      utils::setTxtProgressBar(pb, k);
      k = k + 1

      # set random seed for a given run
      set.seed(random_seeds[j])

      # optimize the model parameters using
      # GenSA algorithm
      par = optimize_parameters(
        par = NULL,
        data = data,
        model = models[i],
        method = method,
        lower = as.numeric(d[1,]),
        upper = as.numeric(d[2,]),
        control = control
      )

      # put optimial parameters in the output
      # matrix
      tmp_parameters[j,1:length(par$par)] = par$par

      # add model output using the estiamted parameters
      predicted_values[j,] = estimate_phenology(
        data = data,
        model = models[i],
        par = par$par
      )
    }

    # stuff everything in a list
    tmp = list("parameters" = tmp_parameters,
               "predicted_values" = predicted_values)

    # append to the list
    model_estimates[i] = list(tmp)
  }

  # close the progress bar element
  close(pb)

  # rename model output for clarity
  names(model_estimates) = models

  # concat comparison data
  comparison = list("modelled" = model_estimates,
                    "measured" = data$transition_dates)

  # check if it's a comparison between two models
  # if so plot the arrrow graph
  if (length(models)==2){

  }

  # return data
  return(comparison)
}

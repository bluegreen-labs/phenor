#' Model validation routine to faciliate model development
#' and quick validation of a (new) model.
#'
#' @param par_ranges: a vector of starting parameter values (function specific)
#' defaults to the parameter ranges as provided with the current models
#' and set forth by Basler (2016)
#' @param model: the model name to be used in optimizing the model
#' @keywords phenology, model, validation
#' @export
#' @examples
#'
#' \dontrun{
#' model_validation(model,par_ranges = "parameter_ranges.csv")
#'
#' # estimate will return the best estimated parameter set given the
#' # validation data
#' }

model_validation = function(model = "TT",
                            dataset = "phenocam_DB",
                            control = list(max.call = 2000),
                            par_ranges = sprintf("%s/extdata/parameter_ranges.csv",
                                                 path.package("phenor"))){

  # if the dataset does not exist
  # in the workspace assume it to be loaded
  # from package storage
  if (is.character(dataset)){
    data(list = dataset)
    dataset = get(dataset)
  }

  # convert to a flat format for speed
  data = flat_format(dataset)

  # read in parameter ranges
  par_ranges = read.table(par_ranges,
                          header = TRUE,
                          sep = ",")

  # subset the parameter range
  if (!any(par_ranges$model == model)){
    stop("parameters are not specified in the standard ")
  }

  # extract parameter ranges is the model is available
  # in the file provided
  d = par_ranges[par_ranges$model == model,]
  d = d[,!is.na(d[1,])]
  d = d[,3:ncol(d)]
  d = as.matrix(d)

  # optimize paramters
  optim.par = optimize_parameters(
    par = NULL,
    data = data,
    model = model,
    method = "GenSA",
    lower = as.numeric(d[1,]),
    upper = as.numeric(d[2,]),
    control = control
  )

  # plot the model output
  val = data$transition_dates
  out = estimate_phenology(
    data = data,
    model = model,
    par = optim.par$par
  )

  RMSE_NULL = sqrt(mean((val - null(data)) ^ 2, na.rm = T))

  RMSE = rmse(par = optim.par$par, data = data, model = model)
  Ac = AICc(measured = val, predicted = out, k = length(optim.par$par))

  plot(val,out,
       main = paste(model,", itterations: ", control$max.call, sep=""),
       xlab = "onset DOY Measured",
       ylab = "onset DOY Modelled",
       pch = 19)
  abline(0,1)
  legend("topleft",legend = sprintf("RMSE: %s",round(RMSE)),bty='n')
  legend("top",legend = sprintf("RMSE NULL: %s",round(RMSE_NULL)),bty='n')
  legend("bottomright",legend = sprintf("AICc: %s",round(Ac$AICc)),bty='n')
  print(summary(lm(val ~ out)))
  return(optim.par$par)
}

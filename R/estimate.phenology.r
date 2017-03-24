#' Calculates the estimated phenophases for a given set of parameters
#' and a specified model (be sure to match parameter and requirements
#' with the model.
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
#' @keywords phenology, model, sequential
#' @export
#' @examples
#' estimate <- estimate.phenology(par,data,model)
#'
#' # estimate will be an estimated timing of a phenophase

# calls a particular model and executes it using a given
# set of parameters and data. Make sure to check the data
# and parameter constraints. Data matrices must match.
estimate.phenology = function(par, data, model, hpc=FALSE){

  # grab the data from the data structure
  temperature = data[[2]]

  # check the size of the data matrix served
  # if the data is a matrix loop over all the
  # rows and calculate the statistics as such
  # if it's just a vector only run once (obviously)
  if(is.matrix(data[[2]])){
    nr_rows = dim(data[[2]])[1]

    # load HPC libraries
    if ( hpc == TRUE ){

      # load HPC libraries, only advisable for very
      # large datasets where unloading and loading of
      # the cluster units does not outweight the speedup
      require(foreach, quietly = TRUE)
      require(doParallel, quietly = TRUE)
      require(parallel, quietly = TRUE)

      numCores <- detectCores()
      cl <- makeCluster(numCores)
      registerDoParallel(cl)

      # compute everything in parallel, make sure your dataset is large
      # enough for this (raster data?)
      results = foreach(i=1:nr_rows,
                        .combine = rbind,
                        .export = model) %dopar% {
        # create a subset of the data
        subset = list(data[[1]],
                      data[[2]][i,],
                      data[[3]][i,],
                      data[[4]][i])
        # call the model, using the model
        do.call(model,list(data=subset,par=par))
      }

      # stop cluster
      stopCluster(cl)

    } else {

      # create output matrix based upon the number
      # of rows in the matrix
      results = matrix(NA,nr_rows,1)

      for ( i in 1:nr_rows ){
        # create a subset of the data
        subset = list(data[[1]],
                      data[[2]][i,],
                      data[[3]][i,],
                      data[[4]][i])

        # call the model, using the model
        results[i,] = do.call(model,list(data=subset,par=par))
      }
    }

  } else {
    results = do.call(model,list(data=subset,par=par))
  }

  # return results
  return(results)
}

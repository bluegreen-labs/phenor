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
#' @param path: the path to daymet tile data formated using the
#' format_daymet() function. All tiles in this directory will be
#' processed with the specified model and data mosaicked and reprojected
#' if requested
#' @param reproject: TRUE / FALSE, reproject to lat-lon (default = TRUE)
#' @keywords phenology, model, post-processing
#' @export
#' @examples
#'
#' \dontrun{
#' # estimate will be an estimated timing of a phenophase
#' estimate = estimate.phenology(par,data,model)
#'
#' # if a path is specified and no other data all tiled
#' # data in this directory will be processed into a map
#' # if a dataset is provided at the data parameter side
#' # only this particular tile will be processed.
#'
#' }

# calls a particular model and executes it using a given
# set of parameters and data. Make sure to check the data
# and parameter constraints. Data matrices must match.

estimate_phenology = function(par,
                              data = NULL,
                              model = "TT",
                              path =  ".",
                              reproject = TRUE) {

  # if the data parameter is not null
  # run the standard procedure
  if(!is.null(data)){

    # convert to a flat format for speed
    data = flat_format(data = data)

    # return model output
    do.call(model, list(data = data, par = par))

  } else {

    # list all rda files
    files = list.files(path = path,
                       pattern = "^phenor_.*\\.rds$",
                       full.names = TRUE)

    # cycle over all files and run the model
    for (i in 1:length(files)){

      # load the rds file (name data)
      data = readRDS(files[i])

      # run the model
      tmp = do.call("TT",list(data = data, par = par))

      # clear the data file
      rm(data)

      # create a phenology map mosaicking all files
      # on the go into one big final map
      if (i == 1){
        phenor_map = tmp
      } else {
        phenor_map = mosaic(phenor_map, tmp, fun = mean)
      }
    }

    # normal output will be lambert conformal conical,
    # if a lat-lon map is required reproject the results
    if (reproject){
      phenor_map = trim(projectRaster(phenor_map, crs = CRS("+init=epsg:4326")))
    }

    # return the map
    return(phenor_map)
  }
}

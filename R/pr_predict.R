#' Predicts phenophases
#'
#' Uses a given set of parameters
#' and a specified model (be sure to match parameter and requirements
#' with the model. Wrapper around a do.call() call for convenience.
#'
#' @param data a nested list of data formatted using one of the format_*()
#' functions
#' @param par a vector of parameter values, this is functions specific
#' @param model the model name to be used in optimizing the model
#' @param path the path to daymet tile data formated using the
#' format_daymet() function. All tiles in this directory will be
#' processed with the specified model and data mosaicked and reprojected
#' if requested
#' @param reproject TRUE / FALSE, reproject to lat-lon (default = TRUE)
#' @param ... extra arguments to pass to the function
#' @return returns point based or gridded estimates of a particular phenopase
#' @keywords phenology, model, post-processing
#' @export
#' @examples
#'
#' \dontrun{
#' # estimate will be an estimated timing of a phenophase
#' estimate <- pr_predict(par, data, model)
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

# pr_predict
pr_predict <- function(
  par,
  data,
  model = "TT",
  path =  ".",
  reproject = TRUE,
  ...
  ){

  # if the data parameter is not null
  # run the standard procedure
  if(!missing(data)){

    # when provided with a file name,
    # load the rds file (name data)
    if(is.character(data)){
      data = readRDS(data)
    }

    # convert to a flat format for speed
    data <- pr_flatten(data = data)

    # return model output
    do.call(model, list(data = data,
                        par = par,
                        ... ))

  } else {

    # list all daymet rds files
    files <- list.files(path = path,
                       pattern = "^phenor_daymet_.*\\.rds$",
                       full.names = TRUE)

    # trap error when no files are found
    if(length(files)==0){
      stop("no Daymet files found for processing...")
    }

    # cycle over all files and run the model
    for (i in 1:length(files)){

      # load the rds file (name data)
      data <- readRDS(files[i])

      # run the model
      tmp <- do.call(model,
                    list(data = data,
                         par = par,
                         ... ))

      # clear the data file
      rm(data)

      # create a phenology map mosaicking all files
      # on the go into one big final map
      if (i == 1){
        phenor_map <- tmp
      } else {
        phenor_map <- terra::mosaic(phenor_map, tmp, fun = mean)
      }
    }

    # normal output will be lambert conformal conical,
    # if a lat-lon map is required reproject the results
    if (reproject){
      phenor_map <-
        terra::project(
          phenor_map,
          "EPSG:4326"
         )
    }

    # return the map
    return(phenor_map)
  }
}

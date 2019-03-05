#' Download Berkeley Earth Gridded mean daily temperature data
#'
#' Global mean daily temperature data.
#'
#' @param path a path where to save the gridded data
#' @param year year to process (requires year - 1 to be present)
#' @return nothing is returned to the R working environment, files are
#' downloaded and stored on disk
#' @keywords phenology, model, data
#' @export
#' @examples
#'
#' # donwload all gridded data for year 2011
#' \dontrun{
#' pr_dl_be(year = 2011)
#'}

# create subset of layers to calculate phenology model output on
pr_dl_be <- function(
  path = tempdir(),
  year = 2011
  ){

  # set server
  server = "http://berkeleyearth.lbl.gov/auto/Global/Gridded"

  # set the decadal splits
  split = seq(1880,as.numeric(format(Sys.Date(),"%Y")),10)

  # download 2 files if on split decade
  if (year %in% split){
    # create download string
    filenames = c(sprintf("Complete_TAVG_Daily_LatLong1_%s.nc",
                          split[max(which(split <= year),na.rm = TRUE) - 1]),
                  sprintf("Complete_TAVG_Daily_LatLong1_%s.nc",
                          split[max(which(split <= year),na.rm = TRUE)]))

  } else {
    # create download string
    filenames = c(sprintf("Complete_TAVG_Daily_LatLong1_%s.nc",
                        split[max(which(split <= year),na.rm = TRUE)]))
  }

  # download missing data
  lapply(filenames, function(i){
    # set download / filename strings
    file_location = sprintf("%s/%s",path,i)
    http_location = sprintf("%s/%s", server, i)

    # try to download the data
    if(!file.exists(file_location)){
      error = try(httr::GET(url = http_location,
                            httr::write_disk(path = file_location,
                                             overwrite = TRUE),
                            httr::progress()),
                  silent = TRUE)
      if (inherits(error, "try-error")){
        file.remove(file_location)
        stop("failed to download the requested data, check your connection")
      }
    } else {
      cat("local file exists, skipping download \n")
    }
  })

  # feedback
  cat("Download complete! \n")
}

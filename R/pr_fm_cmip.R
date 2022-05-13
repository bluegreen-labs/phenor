#' Formatting CMIP6 data
#'
#' Formats CMIP6 data as downloaded through the Copernicus CDS service,
#' using pr_dl_cmip().
#'
#' @param path a path to the gridded data
#' @param year year to process (requires year - 1 to be present / downloaded)
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param model which CMIP6 model to use
#' @param scenario "rcp85", "rcp45", "historical" here rcp covers 2006 - 2100
#' while historical data covers 1950 - 2005
#' @param extent vector with coordinates defining the region of interest defined
#' as xmin, xmax, ymin, ymax in lat/lon (default = c(-80, -70, 40, 50)), note
#' that this differs from the download routine!!
#' @param internal TRUE / FALSE, write data structure to file as RDS
#' (default = FALSE)
#' @return data structure formatted for phenor model optimization and validation
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#' # run with default settings
#' # look for alternative models on the CMIP6
#' # downscaled data page
#' \dontrun{
#' cmip5_data <- pr_fm_cmip()
#'}

# create subset of layers to calculate phenology model output on
pr_fm_cmip <- function(
  path = tempdir(),
  year = 2016,
  offset = 264,
  model = "MIROC6",
  scenario = "ssp5_8_5",
  extent = c(41, -79, 50, -69),
  internal = TRUE
  ) {

  # list all netCDF files in the path
  files <- list.files(path = path,
                     pattern = "*\\.nc",
                     full.names = FALSE)

  # measurements to include in the processing routine
  measurements <- c("tasmax", "tasmin", "pr")

  # download or read data
  data <- lapply(measurements, function(x){

      # filename
      filename <- files[(grepl(x,files) &
                       grepl(gsub("_","",scenario),files) &
                       grepl(toupper(model),toupper(files))
                       )]

      # if the file exist use the local file
      if (length(filename) != 0){
        r <- terra::rast(file.path(path, filename), x)

        # extract dates from meta-data
        dates <- terra::time(r)

        # define start and end date
        start <- as.POSIXct(sprintf("%s-01-01", year-1), tz = "UTC")
        end <- as.POSIXct(sprintf("%s-12-31", year), tz = "UTC")

        # find locations of the layers to select
        layers <- which(dates >= start & dates <= end)

        # only return valid layers
        r <- terra::subset(r, layers)
        return(r)
      } else {
        stop("Required files not available: Check your path variable!")
      }
   })

  # crop the data if large files are provided this saves
  # time and data loads

  # set order of extent element for terra
  e <- c(
    extent[2],
    extent[4],
    extent[3],
    extent[1]
  )

  Tmini <- terra::crop(data[[2]], terra::ext(e)) - 273.15
  Tmaxi <- terra::crop(data[[1]], terra::ext(e)) - 273.15
  Pi <- terra::crop(data[[3]], terra::ext(e))
  temp <- (Tmaxi + Tmini) / 2

  # extract the yday and year strings
  # year strings, depends on how things are subset
  # and pasted back together
  dates <- terra::time(Tmini)
  yday <- as.numeric(format(dates,"%j"))
  years <- as.numeric(format(dates,"%Y")) - year + 2

  # calculate if the previous year was a leap year
  # to account for this offset
  leap_year = ifelse((year-1%%4==0 & year-1%%100!=0) | year-1%%400==0,
                     TRUE,
                     FALSE)

  # select layers to subset using this year and yday data
  # account for leap years included in the NEX data
  if(leap_year){
    layers = which((years == 1 & yday >= offset) |
                   (years == 2 & yday < (offset - 1)))
  } else {
    layers = which((years == 1 & yday >= offset) |
                     (years == 2 & yday < offset))
  }

  # check if all layers are available in the dataset
  if (length(layers) < 365){
   stop("The selected dataset does not cover your data range!")
  }

  # subset raster data based upon the offset chosen
  temp <- terra::subset(temp, layers)
  Tmini <- terra::subset(Tmini, layers)
  Tmaxi <- terra::subset(Tmaxi, layers)
  Pi <- terra::subset(Pi, layers)

  # extract georeferencing info to be passed along
  ext <- terra::ext(temp)
  proj <- terra::crs(temp, proj = TRUE)
  size <- dim(temp)

  # grab coordinates
  location <- terra::xyFromCell(temp, 1:prod(dim(temp)[1:2]))

  # create doy vector
  if (offset < length(layers)){
    doy = c(offset:length(layers),1:(offset - 1))
  } else {
    doy = 1:length(layers)
  }

  # create daylength matrix
  Li <- lapply(location[,2],
              FUN = function(x){
                  daylength(doy = doy, latitude = x)
                })
  Li <- t(do.call("rbind", Li))

  # recreate the validation data structure (new format)
  # but with concatted data
  data = list("site" = NULL,
              "location" = location,
              "doy" = doy,
              "transition_dates" = NULL,
              "Ti" = t(terra::as.matrix(temp)),
              "Tmini" = t(terra::as.matrix(Tmini)),
              "Tmaxi" = t(terra::as.matrix(Tmaxi)),
              "Li" = Li,
              "Pi" = t(terra::as.matrix(Pi)),
              "VPDi" = NULL,
              "georeferencing" = list("extent" = as.vector(ext),
                                      "projection" = proj,
                                      "size" = size)
  )

  # assign a class for post-processing
  class(data) = "phenor_map_data"

  # return the formatted, faster data format
  # either internally or saved as an rds (binary R data file)
  if (internal){
    return(data)
  } else {
    saveRDS(data, file = sprintf("%s/phenor_cmip_data_%s_%s_%s.rds",
                                 path,
                                 model,
                                 year,
                                 scenario))
    # clean up
    gc()
  }
}

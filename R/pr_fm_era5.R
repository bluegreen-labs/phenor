#' Formatting ERA5 data
#'
#' Formats ERA5 data as downloaded through the Copernicus CDS service,
#' using pr_dl_era5().
#'
#' @param path a path to the gridded data
#' @param file a netCDF containing ERA5 data, by default this will be a
#'  file named era5.nc. Alternatively you can point to a different file
#'  in the path (e.g. when you downloaded multiple source files)
#' @param year year to process (requires year - 1 to be present / downloaded)
#' @param offset offset of the time series in DOY (default = 264, sept 21)
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
pr_fm_era5 <- function(
  path = tempdir(),
  file = "era5.nc",
  year = 2016,
  offset = 264,
  internal = TRUE,
  extent
  ) {

  # measurements to include in the processing routine
  # might increase the number of variables later to
  # allow for more flexible models, however most models
  # are only driven by temperature and externally calculated
  # daylength values
  measurements <- c("t2m", "tp")

  # missing doesn't work in nested
  # function calls, create own variable
  crop_nc <- !missing(extent)

  # download or read data
  data <- lapply(measurements, function(x){

      # if the file exist use the local file
      if (file.exists(file.path(path, file))){

        # load in the raster file
        r <- terra::rast(file.path(path, file), x)

        # crop data to size
        if(crop_nc){
          r <- try(terra::crop(r, extent))

          # trap cropping errors
          if(inherits(r, "try-error")){
            stop("Failed cropping to your extent, please check coordinate
                 order and ranges")
          }
        }

        # extract dates from meta-data
        # and reformat to a year-doy format
        # for aggregation on a daily level
        dates <- terra::time(r)
        year_doy <- format(dates, "%Y-%j")

        # define start and end date
        start <- as.POSIXct(sprintf("%s-01-01", year-1), tz = "UTC")
        end <- as.POSIXct(sprintf("%s-12-31", year), tz = "UTC")

        # find locations of the layers to select
        layers <- which(dates >= start & dates <= end)

        # subset the layers
        # only return valid layers
        r <- terra::subset(r, layers)

        # aggregate to a daily level for either
        # temperature or precipitation values
        if(x == "t2m"){

          # aggregate by index year_doy
          Tmaxi <- terra::tapp(r, index = year_doy, fun = max)
          Tmini <- terra::tapp(r, index = year_doy, fun = min)
          temp <- terra::tapp(r, index = year_doy, fun = mean)

          # add in the time stamps
          terra::time(Tmaxi) <- unique(as.Date(dates))
          terra::time(Tmini) <- unique(as.Date(dates))
          terra::time(temp) <- unique(as.Date(dates))

          return(list(Tmaxi, Tmini, temp))

        } else {
          Pi <- terra::tapp(r, index = year_doy, fun = sum)
          terra::time(Pi) <- unique(as.Date(dates))
          return(Pi)
        }

      } else {
        stop("Required files not available: Check your path variable!")
      }
   })

  # split out data from nested
  # structure
  Tmini <- data[[1]][[2]]
  Tmaxi <- data[[1]][[1]]
  Pi <- data[[2]]
  temp <- data[[1]][[3]]

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
  # note this can be moved up into the processing
  # above (more efficient probably!!)
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
              "Ti" = t(terra::as.matrix(temp) - 273.15),
              "Tmini" = t(terra::as.matrix(Tmini) - 273.15),
              "Tmaxi" = t(terra::as.matrix(Tmaxi) - 273.15),
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

#' Preprocessing of NASA Earth Exchange Global Daily Downscaled Projections
#' (NEX-GDDP).
#'
#' @param path a path to the gridded data (original max / min temp. / precip files)
#' @param year year to process (requires year - 1 to be present / downloaded)
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param model which CMIP5 model to use
#' @param scenario "rcp85", "rcp45", "historical" here rcp covers 2006 - 2100
#' while historical data covers 1950 - 2005
#' @param extent vector with coordinates defining the region of interest defined
#' as xmin, xmax, ymin, ymax in lat/lon (default = c(-74,-65, 40, 48))
#' @param internal TRUE / FALSE, write data structure to file as RDS
#' (default = FALSE)
#' @return data structure formatted for phenor model optimization and validation
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#' # run with default settings
#' # look for alternative models on the CMIP5
#' # downscaled data page
#' \dontrun{
#' cmip5_data = format_cmip5()
#'}

# create subset of layers to calculate phenology model output on
format_cmip5 = function(path = "~",
                        year = 2016,
                        offset = 264,
                        model = "IPSL-CM5A-MR",
                        scenario = "rcp85",
                        extent = c(-128,-65, 24, 50),
                        internal = TRUE){

  # list all netCDF files in the path
  files = list.files(path = path,
                     pattern = "*\\.nc",
                     full.names = FALSE)

  # measurements to include in the processing routine
  measurements = c("tasmax", "tasmin", "pr")

  # loop over the two years needed
  # depending on the offset
  for (i in c(year-1, year)){

    # download or read data
    data = lapply(measurements, function(x){

      # filename
      filename = files[which(grepl(i,files) &
                       grepl(x,files) &
                       grepl(scenario,files) &
                       grepl(toupper(model),files)
                       )]

      # if the file exist use the local file
      if (length(filename) != 0){
        r = raster::brick(sprintf("%s/%s", path, filename))
        return(r)
      } else {
        stop("Required files not available!
  Data will not be downloaded automatically due to large file sizes.
  Please download data using the download_cmip5() function.")
      }
    })

    # stack the temperature data to take the mean
    # using stackApply()
    temp_data_stack = raster::stack(data[[1]],
                              data[[2]])

    # shift the data, longitudes run from 0 - 360
    temp_data_stack = raster::shift(temp_data_stack, x = -360)
    min_temp = raster::shift(data[[2]], x = -360)
    max_temp = raster::shift(data[[1]], x = -360)
    precip = raster::shift(data[[3]], x = -360)

    # crop data for faster processing
    temp_data_stack = raster::crop(temp_data_stack, raster::extent(extent))
    min_temp = raster::crop(min_temp, raster::extent(extent))
    max_temp = raster::crop(max_temp, raster::extent(extent))
    precip = raster::crop(precip, raster::extent(extent))

    # calculate mean temperature on cropped data for speed
    l = raster::nlayers(temp_data_stack)/2
    mean_temp = raster::stackApply(temp_data_stack,
                                   indices = rep(1:l,2),
                                   fun = mean,
                                   na.rm = TRUE)

    # combine data in one big stack
    if (year != i){
        temp = mean_temp
        Tmini = min_temp
        Tmaxi = max_temp
        Pi = precip
    } else {
      # combine data with previous year's
      temp = raster::stack(temp, mean_temp)
      Tmini = raster::stack(Tmini, min_temp)
      Tmaxi = raster::stack(Tmaxi, max_temp)
      Pi = raster::stack(Pi, precip)
    }

    # clear intermediate variables
    rm(list = c("mean_temp","min_temp","max_temp","precip"))
  }

  # extract the yday and year strings
  dates = as.Date(names(Tmini),"X%Y.%m.%d")
  yday = as.numeric(format(dates,"%j"))
  years = as.numeric(format(dates,"%Y"))

  # calculate if the previous year was a leap year
  # to account for this offset
  leap_year = ifelse((year-1%%4==0 & year-1%%100!=0) | year-1%%400==0,
                     TRUE,
                     FALSE)

  # select layers to subset using this year and yday data
  # account for leap years included in the NEX data
  if(leap_year){
    layers = which((years == (year - 1) & yday >= offset) |
                   (years == year & yday < (offset - 1)))
  } else {
    layers = which((years == (year - 1) & yday >= offset) |
                     (years == year & yday < offset))
  }

  # check if all layers are available in the dataset
  if (length(layers) < 365){
   stop("The selected dataset does not cover your data range!")
  }

  # subset raster data
  temp = raster::subset(temp, layers)
  Tmini = raster::subset(Tmini, layers)
  Tmaxi = raster::subset(Tmaxi, layers)
  Pi = raster::subset(Pi, layers)

  # extract georeferencing info to be passed along
  ext = raster::extent(temp)
  proj = raster::projection(temp)
  size = dim(temp)

  # grab coordinates
  location = sp::SpatialPoints(
    sp::coordinates(temp),
    proj4string = sp::CRS(proj))
  location = t(sp::spTransform(location,
                               sp::CRS("+init=epsg:4326"))@coords[,2:1])

  # create doy vector
  if (offset < length(layers)){
    doy = c(offset:length(layers),1:(offset - 1))
  } else {
    doy = 1:length(layers)
  }

  # create daylength matrix
  Li = lapply(location[1,],
              FUN = function(x){
                  daylength(doy = doy, latitude = x)
                })
  Li = t(do.call("rbind", Li))

  # recreate the validation data structure (new format)
  # but with concatted data
  data = list("site" = NULL,
              "location" = location,
              "doy" = doy,
              "transition_dates" = NULL,
              "Ti" = t(raster::as.matrix(temp)) - 273.15,
              "Tmini" = t(raster::as.matrix(Tmini)) - 273.15,
              "Tmaxi" = t(raster::as.matrix(Tmaxi)) - 273.15,
              "Li" = Li,
              "Pi" = t(raster::as.matrix(Pi)),
              "VPDi" = NULL,
              "georeferencing" = list("extent" = ext,
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
    saveRDS(data, file = sprintf("%s/phenor_cmip5_data_%s_%s_%s.rds",
                                 path,
                                 model,
                                 year,
                                 scenario))
    # clean up
    gc()
  }
}

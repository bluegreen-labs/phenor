#' Preprocessing CMIP5 model runs as hosted by:
#' http://gdo-dcp.ucllnl.org/downscaled_cmip_projections/dcpInterface.html
#'
#' The function will format the data into the correct phenor structure for
#' post-processing.
#'
#' NOTE:
#' Although one can access the 1/16th degree data, by default the 1x1 degree
#' data should be preferred (due to preprocessing load). Furthermore, certain
#' firewalls will block timely access to the data which will result in failure
#' to proceed. In this case download the data first from the website and process
#' it by directing the function to the correct download path.
#'
#' @param path a path to the gridded data (original max / min temp. files)
#' @param year year to process (requires year - 1 to be present / downloaded)
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param model which CMIP5 model to use
#' @param resolution "1x1" or "16th", "16th" will result in a significant
#' calculation overhead
#' @param scenario "rcp85", "rcp45", "historical" here rcp covers 2006 - 2100
#' while historical data covers 1950 - 2005
#' @param internal TRUE / FALSE, write data structure to file or not (as .rds)
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
                        resolution = "1x1",
                        scenario = "rcp85",
                        internal = TRUE){
  # set server
  server = "ftp://gdo-dcp.ucllnl.org/pub/dcp/archive/cmip5/loca/LOCA_2016-04-02"

  # measurements to average
  measurements = c("tasmax","tasmin","pr")

  # loop over the two years needed
  # depending on the offset
  for (i in c(year-1, year)){

    # download or read data
    data = lapply(measurements, function(x){

      # filename
      filename = sprintf("%s_day_%s_%s_r1i1p1_%s0101-%s1231.LOCA_2016-04-02.%s.nc",
                         x,
                         model,
                         scenario,
                         i,
                         i,
                         resolution)

      # if the file exist use the local file
      if (file.exists(sprintf("%s/%s", path, filename))){
        r = raster::brick(sprintf("%s/%s", path, filename))
        return(r)
      } else {
        # url
        url = sprintf("%s/%s/%s/%s/r1i1p1/%s/%s",
                      server,
                      model,
                      resolution,
                      scenario,
                      x,
                      filename)

        # download and return the data
        error = try(curl::curl_download(url = url,
                                        destfile = sprintf("%s/%s",
                                                           path,
                                                           filename),
                                        quiet = TRUE))

        if (inherits(error, "try-error")){
          file.remove(sprintf("%s/%s", path, filename))
          stop('Server not reachable, remember R does not support passive FTP!')
        } else {
          r = raster::brick(sprintf("%s/%s", path, filename))
        }
        return(r)
      }
    })

    # stack the temperature data to take the mean
    # using stackApply()
    temp_data_stack = raster::stack(data[[1]],
                              data[[2]])

    l = raster::nlayers(temp_data_stack)/2
    mean_temp = raster::stackApply(temp_data_stack,
                                   indices = rep(1:l,2),
                                   fun = mean,
                                   na.rm = TRUE)

    # shift the data, longitudes run from 0 - 360
    mean_temp = raster::shift(mean_temp, x = -360)
    min_temp = raster::shift(data[[2]], x = -360)
    max_temp = raster::shift(data[[1]], x = -360)
    precip = raster::shift(data[[3]], x = -360)

    # drop layer 366 on leap year
    if (raster::nlayers(mean_temp) == 366){
      mean_temp = raster::dropLayer(mean_temp, 366)
      min_temp = raster::dropLayer(min_temp, 366)
      max_temp = raster::dropLayer(max_temp, 366)
      precip = raster::dropLayer(precip, 366)
    }

    # combine data in one big stack
    if (year != i){
      temp = mean_temp
      Tmini = min_temp
      Tmaxi = max_temp
      Pi = precip
    } else {
      temp = raster::stack(temp, mean_temp)
      Tmini = raster::stack(Tmini, min_temp)
      Tmaxi = raster::stack(Tmaxi, max_temp)
      Pi = raster::stack(Pi, precip)
    }

    # clear variable
    rm(list = c("mean_temp","min_temp","max_temp","precip"))
  }

  # extract the yday and year strings
  yday = floor(as.numeric(unlist(lapply(strsplit(names(temp),"_"),"[[",2))))
  years = unlist(lapply(strsplit(names(temp),"\\."),"[[",2))

  # select layers to subset using this year and yday data
  layers = which((years == 1 & yday >= offset) |
                   (years == 2 & yday < offset))

  # check if all layers are available in the dataset
  # needs a fix for crossing boundaries between decadal subsets !!
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
  # either internally or saved as an rda (binary R data file)
  if (internal){
    return(data)
  } else {
    saveRDS(data, file = sprintf("%s/phenor_cmip5_data_%s.rds",path, year))
  }
}

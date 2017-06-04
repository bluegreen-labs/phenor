#' Preprocessing CMIP5 model runs as hosted by:
#' http://gdo-dcp.ucllnl.org/downscaled_cmip_projections/dcpInterface.html
#'
#' Although one can access the 1/16th degree data, by default the 1x1 degree
#' data should be preferred (due to preprocessing load).
#'
#' @param path: a path to the gridded data (original max / min temp. files)
#' @param year: year to process (requires year - 1 to be present / downloaded)
#' @param offset: offset of the time series in DOY (default = 264, sept 21)
#' @param model: which CMIP5 model to use
#' @param resolution: "1x1" or "16th", "16th" will result in a significant
#' calculation overhead
#' @param scenario: "rcp85", "rcp45", "historical" here rcp covers 2006 - 2100
#' while historical data covers 1950 - 2005
#' @param internal: TRUE / FALSE, write data structure to file or not (as .rds)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples

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
  measurements = c("tasmax","tasmin")

  # loop over the two years needed
  # depending on the offset
  for (i in c(year-1, year)){

    # download or read data
    temp_data = lapply(measurements,function(x){

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
                                        destfile = sprintf("%s/%s", path, filename)))

        if (inherits(error, "try-error")){
          file.remove(sprintf("%s/%s", path, filename))
          stop('File not found')
        } else {
          r = raster::brick(sprintf("%s/%s", path, filename))
        }
        return(r)
      }
    })

    # stack the temperature data to take the mean
    # using stackApply()
    temp_data_stack = raster::stack(temp_data[[1]],
                              temp_data[[2]]) - 273.15
    assign("test",temp_data_stack, envir = .GlobalEnv)

    l = nlayers(temp_data_stack)/2
    mean_temp = raster::stackApply(temp_data_stack,
                                   indices = rep(1:l,2),
                                   fun = mean,
                                   na.rm = TRUE)

    # shift the data, longitudes run from 0 - 360
    mean_temp = raster::shift(mean_temp, x = -360)

    # drop layer 366 on leap year
    if (nlayers(mean_temp) == 366){
      mean_temp = raster::dropLayer(mean_temp, 366)
    }

    # combine data in one big stack
    if (year != i){
      temperature = mean_temp
    } else {
      temperature = stack(temperature, mean_temp)
    }

    # clear variable
    rm(mean_temp)
  }

  # extract the yday and year strings
  yday = floor(as.numeric(unlist(lapply(strsplit(names(temperature),"_"),"[[",2))))
  years = unlist(lapply(strsplit(names(temperature),"\\."),"[[",2))

  # select layers to subset using this year and yday data
  layers = which((years == 1 & yday >= offset) |
                   (years == 2 & yday < offset))

  # check if all layers are available in the dataset
  # needs a fix for crossing boundaries between decadal subsets !!
  if (length(layers) < 365){
    stop("The selected dataset does not cover your data range!")
  }

  # subset raster data, correct for Kelvin scale
  temperature = subset(temperature, layers)

  # shift data when offset is < 365
  if (offset < length(layers)){
    doy = c(offset:length(layers),1:(offset - 1))
  } else {
    doy = 1:length(layers)
  }

  # convert temperature data to matrix
  Ti = t(as.matrix(temperature))

  # extract georeferencing info to be passed along
  ext = extent(temperature)
  proj = projection(temperature)
  size = dim(temperature)

  # grab coordinates
  location = SpatialPoints(coordinates(temperature),
                           proj4string = CRS(proj))
  location = t(spTransform(location, CRS("+init=epsg:4326"))@coords[,2:1])

  # create daylength matrix
  Li = matrix(rep(1:length(layers), prod(size[1:2])),
              length(layers),
              prod(size[1:2]))

  latitude = matrix(location[1,],
                    length(layers),
                    prod(size[1:2]), byrow = TRUE)

  # calculate daylength
  Li = daylength(doy = Li, latitude = latitude)[1]
  Li = Li[[1]]

  # create doy vector
  if (offset < length(layers)){
    doy = c(offset:length(layers),1:(offset - 1))
  } else {
    doy = 1:length(layers)
  }

  # recreate the validation data structure (new format)
  # but with concatted data
  data = list("site" = NULL,
              "location" = location,
              "doy" = doy,
              "transition_dates" = NULL,
              "Ti" = Ti,
              "Li" = Li,
              "Pi" = NULL,
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
    saveRDS(data, file = sprintf("%s/phenor_data_%s_%s.rds",path, year, tile))
  }
}

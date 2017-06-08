#' Preprocessing of Berkeley Earth Gridded temperature data
#'
#' @param path: a path to the gridded data, only average temperature
#' data will be considered
#' @param year: year to process (requires year - 1 to be present)
#' @param offset: offset of the time series in DOY (default = 264, sept 21)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples

# create subset of layers to calculate phenology model output on
format_berkeley_earth = function(path = "~",
                                 year = 2011,
                                 offset = 264,
                                 bounding_box = c(-126, -66, 23, 54),
                                 internal = TRUE){

  # set server
  server = "http://berkeleyearth.lbl.gov/auto/Global/Gridded"

  # set the decadal splits
  split = seq(1880,as.numeric(format(Sys.Date(),"%Y")),10)

  # download 2 files if on split decade
  if (year %in% split){
    # create download string
    filenames = c(sprintf("Complete_TAVG_Daily_LatLong1_%s.nc",
                          split[max(which(split <= year),na.rm=TRUE) - 1]),
                  sprintf("Complete_TAVG_Daily_LatLong1_%s.nc",
                          split[max(which(split <= year),na.rm=TRUE)]))

  } else {
    # create download string
    filenames = sprintf("Complete_TAVG_Daily_LatLong1_%s.nc",
                        split[max(which(split <= year),na.rm=TRUE)])
  }

  # download missing data
  for (i in filenames){

    # set download / filename strings
    file_location = sprintf("%s/%s",path,i)
    http_location = sprintf("%s/%s", server, i)

    if(!file.exists(file_location)){
      # try to download the data
      error = try(curl::curl_download(http_location,
                                      file_location,
                                      mode="w",
                                      quiet=TRUE),silent=TRUE)
      if (inherits(error, "try-error")){
        # remove file header
        file.remove(file_location)
        stop("failed to download the requested data, check your connection")
      }
    } else {
      cat("local file exists, skipping download \n")
    }
  }

  # grab raster data from netcdf files
  climatology = raster::brick(sprintf("%s/%s",path,filenames[1]),
                              varname = "climatology")
  delta = do.call("stack",
                  lapply(filenames, FUN = function(x){
                    raster::brick(sprintf("%s/%s",path,x),
                                  varname = "temperature")}))

  # get date elements from netcdf files
  years = do.call("c",lapply(filenames,FUN = function(x){
    nc = ncdf4::nc_open(sprintf("%s/%s",path,x))
    ncdf4::ncvar_get(nc,"year")
  }))
  yday = do.call("c",lapply(filenames,FUN = function(x){
    nc = ncdf4::nc_open(sprintf("%s/%s",path,x))
    ncdf4::ncvar_get(nc,"day_of_year")
  }))

  # select layers to subset
  layers = which((years == (year - 1) & yday >= offset) |
                   (years == year & yday < offset))

  # check if all layers are available in the dataset
  # needs a fix for crossing boundaries between decadal subsets !!
  if (length(layers) < 365){
    stop("The selected Berkeley Earth dataset does not cover your data range!")
  }

  # subset raster data, do this before cropping to make things faster
  delta = raster::subset(delta, layers)

  # crop data if necessary
  if (!is.null(bounding_box)){
    if(length(bounding_box)==4){
      climatology = raster::crop(climatology,extent(bounding_box))
      delta = raster::crop(delta,extent(bounding_box))
    } else {
      stop("not enough coordinate points to properly crop data!")
    }
  }

  # shift data when offset is < 365
  if (offset < length(layers)){
    doy = c(offset:length(layers),1:(offset - 1))
  } else {
    doy = 1:length(layers)
  }

  # calculate absolute temperatures not differences
  # with the mean
  temperature = raster::stackApply(raster::stack(delta,climatology),
                           indices = c(doy,1:length(layers)),
                           fun = sum )
  temperature[temperature == 0] = NA

  # convert temperature data to matrix use the raster as.matrix() function!!
  Ti = t(raster::as.matrix(temperature))

  # extract georeferencing info to be passed along
  ext = raster::extent(temperature)
  proj = raster::projection(temperature)
  size = dim(temperature)

  cat("calculating daylength \n")

  # grab coordinates
  location = sp::SpatialPoints(coordinates(temperature),
                           proj4string = CRS(proj))
  location = t(sp::spTransform(location, CRS("+init=epsg:4326"))@coords[,2:1])

  # create doy vector
  if (offset < length(layers)){
    doy = c(offset:length(layers),1:(offset - 1))
  } else {
    doy = 1:length(layers)
  }

  # create daylength matrix
  Li = lapply(location[1,],
              FUN = function(x){
                unlist(daylength(doy = doy, latitude = x)[1])
              })
  Li = t(do.call("rbind",Li))

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


  # return the formatted, faster data format
  # either internally or saved as an rda (binary R data file)
  if (internal){
    return(data)
  } else {
    saveRDS(data, file = sprintf("%s/phenor_be_data_%s_%s.rds",path, year))
  }
}

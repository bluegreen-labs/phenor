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
format_berkeley_earth = function(file = "~/Downloads/Complete_TAVG_Daily_LatLong1_2010.nc",
                               year = 2011,
                               offset = 264,
                               bounding_box = c(-180,-40,0,90),
                               internal = TRUE){

  # some feedback
  cat("calculating average daily temperatures, or load from file \n")

  # read in netCDF data using the raster package
  climatology = brick(file, varname = "climatology") # baseline averages
  delta = brick(file, varname = "temperature") # anomaly

  # crop data if necessary
  if (!is.null(bounding_box)){
    if(length(bounding_box)==4){
      climatology = crop(climatology,extent(bounding_box))
      delta = crop(delta,extent(bounding_box))
    } else {
      stop("not enough coordinate points to properly crop data!")
    }
  }

  # get dates from ncdf layers
  nc = ncdf4::nc_open(file)
  years = ncdf4::ncvar_get(nc,"year")
  yday = ncdf4::ncvar_get(nc,"day_of_year")

  # select layers to subset
  layers = which((years == (year - 1) & yday >= offset) |
                   (years == year & yday < offset))

  # check if all layers are available in the dataset
  # needs a fix for crossing boundaries between decadal subsets !!
  if (length(layers)!=365){
    stop("The selected Berkeley Earth dataset does not cover your data range!")
  }

  # subset raster data
  delta = subset(delta, layers)

  # shift data when offset is < 365
  if (offset < 365){
    doy = c(offset:365,1:(offset - 1))
  } else {
    doy = 1:365
  }

  # calculate absolute temperatures not differences
  # with the mean
  temperature = stackApply(stack(delta,climatology),
             indices = c(doy,1:365),
             fun = sum )
  temperature[temperature == 0] = NA

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
  Li = matrix(rep(1:365, prod(size[1:2])),
              365,
              prod(size[1:2]))

  latitude = matrix(location[1,],
                    365,
                    prod(size[1:2]), byrow = TRUE)

  # calculate daylength
  Li = daylength(doy = Li, latitude = latitude)[1]
  Li = Li[[1]]

  # create doy vector
  if (offset < 365){
    doy = c(offset:365,1:(offset - 1))
  } else {
    doy = 1:365
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


  # return the formatted, faster data format
  # either internally or saved as an rda (binary R data file)
  if (internal){
    return(data)
  } else {
    saveRDS(data, file = sprintf("%s/phenor_data_%s_%s.rds",path, year, tile))
  }
}

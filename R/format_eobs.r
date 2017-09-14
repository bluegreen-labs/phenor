#' Preprocessing E-OBS data as hosted by:
#' http://www.ecad.eu/
#'
#' Haylock, M.R., N. Hofstra, A.M.G. Klein Tank, E.J. Klok, P.D.
#' Jones, M. New. 2008: A European daily high-resolution gridded dataset
#' of surface temperature and precipitation.
#' J. Geophys. Res (Atmospheres), 113, D20119
#'
#' Please register before downloading data, unzip data into uncompressed
#' netCDF files before processing.
#'
#' @param path: a path of the gridded data netCDF files
#' @param year: year to process (requires year - 1 to be present / downloaded)
#' @param offset: offset of the time series in DOY (default = 264, sept 21)
#' @param resolution: 0.25 or 0.50 degree data, only the normal rectangular data
#' will be processed not the polar variety.
#' calculation overhead
#' @param internal: TRUE / FALSE, write data structure to file or not (as .rds)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' # run with default settings netCDF E-OBS files will be
#' # looked for in your home directory
#' eobs_data = format_eobs()
#'}

# create subset of layers to calculate phenology model output on
format_eobs = function(path = "~/tmp/eobs/",
                        year = 2014,
                        offset = 264,
                        resolution = "0.25",
                        internal = TRUE){

  # measurements to average
  measurements = c("tg","rr")

  # download or read data
  eobs_data = lapply(measurements,function(x){

    # filename
    filename = sprintf("%s_%sdeg_reg_v15.0.nc",
                    x,
                    resolution)

    # if the file exist use the local file
    if (file.exists(sprintf("%s/%s", path, filename))){
      r = raster::brick(sprintf("%s/%s", path, filename))
      return(r)
    } else {
      stop('No E-OBS files found in the referred path !')
    }
  })

  # extract the yday and year strings and convert to numbers
  eobs_date = as.Date(eobs_data[[1]]@z$Date)
  yday = as.numeric(format(as.Date(eobs_data[[1]]@z$Date),"%j"))
  years = as.numeric(format(as.Date(eobs_data[[1]]@z$Date),"%Y"))

  # select layers to subset using this year and yday data
  layers = which((years == (year - 1) & yday >= offset) |
                   (years == year & yday < offset))

  # check if all layers are available in the dataset
  # needs a fix for crossing boundaries between decadal subsets !!
  if (length(layers) < 365){
    stop("The selected dataset does not cover your data range!")
  }

  # subset raster data, correct for Kelvin scale
  temperature = raster::subset(eobs_data[[1]], layers)
  #precipitation = raster::subset(temperature, layers)

  # shift data when offset is < 365
  if (offset < length(layers)){
    doy = c(offset:length(layers),1:(offset - 1))
  } else {
    doy = 1:length(layers)
  }

  # convert temperature data to matrix
  Ti = t(raster::as.matrix(temperature))

  # extract georeferencing info to be passed along
  ext = extent(temperature)
  proj = projection(temperature)
  size = dim(temperature)

  # grab coordinates
  location = SpatialPoints(coordinates(temperature),
                           proj4string = CRS(proj))
  location = t(spTransform(location, CRS("+init=epsg:4326"))@coords[,2:1])

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

  # assign a class for post-processing
  class(data) = "phenor_map_data"

  # return the formatted, faster data format
  # either internally or saved as an rda (binary R data file)
  if (internal){
    return(data)
  } else {
    saveRDS(data, file = sprintf("%s/phenor_eobs_data_%s.rds",path, year))
  }

}

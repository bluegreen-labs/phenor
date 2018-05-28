#' Preprocessing E-OBS data as hosted by:
#' http://www.ecad.eu/
#'
#' Download E-OBS data as described by:
#' Haylock, M.R., N. Hofstra, A.M.G. Klein Tank, E.J. Klok, P.D.
#' Jones, M. New. 2008: A European daily high-resolution gridded dataset
#' of surface temperature and precipitation.
#' J. Geophys. Res (Atmospheres), 113, D20119
#'
#' Please register before downloading data, unzip data into uncompressed
#' netCDF files before processing.
#'
#' @param path a path of the gridded data netCDF files
#' @param year year to process (requires year - 1 to be present / downloaded)
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param resolution 0.25 or 0.50 degree data, only the normal rectangular data
#' will be processed not the polar variety.
#' calculation overhead
#' @param internal TRUE / FALSE, write data structure to file or not (as .rds)
#' @return Returns spatial model input data from long term E-OBS modelled data.
#' This data can be run by models specified in the phenor package.
#' @keywords phenology, model, preprocessing, climate data
#' @export
#' @examples
#'
#' \dontrun{
#' # run with default settings netCDF E-OBS files will be
#' # looked for in your home directory
#' eobs_data = format_eobs()
#'}

format_eobs = function(path = "~",
                        year = 2014,
                        offset = 264,
                        resolution = 0.25,
                        internal = TRUE){

  # download or read data
  eobs_data = lapply( c("tg","rr","elev"),function(x){

    # filename
    filename = sprintf("%s_%sdeg_reg[^/]*\\.nc",
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

  # subset raster data
  Ti = t(raster::as.matrix(raster::subset(eobs_data[[1]], layers)))
  Pi = t(raster::as.matrix(raster::subset(eobs_data[[2]], layers)))
  altitude = t(raster::as.matrix(eobs_data[[3]]$layer))

  # extract georeferencing info to be passed along
  ext = raster::extent(eobs_data[[1]])
  proj = raster::projection(eobs_data[[1]])
  size = dim(eobs_data[[1]])

  # grab coordinates
  location = sp::SpatialPoints(sp::coordinates(eobs_data[[1]]),
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
  Li = t(do.call("rbind",Li))

  # recreate the validation data structure (new format)
  # but with concatted data
  data = list("site" = NULL,
              "location" = location,
              "doy" = doy,
              "transition_dates" = NULL,
              "Ti" = Ti,
              "Tmini" = NULL,
              "Tmaxi" = NULL,
              "Li" = Li,
              "Pi" = Pi,
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
    saveRDS(data, file = sprintf("%s/phenor_eobs_data_%s.rds",path, year))
  }
}

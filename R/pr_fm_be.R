#' Preprocessing of Berkeley Earth Gridded temperature data
#'
#' @param path a path to the gridded data, only average temperature
#' data will be considered
#' @param year year to process (requires year - 1 to be present)
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param extent geographic coordinates constraining the output, defined
#' as bottom left, top right c(lon1, lat1, lon2, lat2) (default =
#' c(-126, -66, 23, 54))
#' @param internal return results as an R variable or save as an RDS file
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' # run with default settings
#' # extracts values of the referenced publication
#' # see github README
#' be_data <- pr_fm_be()
#'}

# create subset of layers to calculate phenology model output on
pr_fm_be <- function(
  path = tempdir(),
  year = 2011,
  offset = 264,
  extent = c(-126,-66, 23, 54),
  internal = TRUE
  ) {

  # download all data if it doesn't exist
  pr_dl_be(path = path, year = year)

  # set the decadal splits
  split <- seq(1880, as.numeric(format(Sys.Date(),"%Y")), 10)

  # create download string
  # download 2 files if on split decade
  if (year %in% split){
    tavg_files <- c(sprintf("Complete_TAVG_Daily_LatLong1_%s.nc",
                           split[max(which(split <= year),na.rm=TRUE) - 1]),
                   sprintf("Complete_TAVG_Daily_LatLong1_%s.nc",
                           split[max(which(split <= year),na.rm=TRUE)]))
  } else {
    filenames <- c(sprintf("Complete_TAVG_Daily_LatLong1_%s.nc",
                          split[max(which(split <= year),na.rm=TRUE)]))
  }

  # grab raster data from netcdf files, with the climatology
  # the long term mean and delta the differences
  climatology <- raster::brick(sprintf("%s/%s",path,filenames[1]),
                              varname = "climatology")
  delta <- do.call("stack",
                  lapply(filenames, FUN = function(x){
                    raster::brick(sprintf("%s/%s",path,x),
                                  varname = "temperature")}))

  # get date elements from netcdf files
  years <- do.call("c",lapply(filenames,FUN = function(x){
    nc <- ncdf4::nc_open(sprintf("%s/%s",path,x))
    ncdf4::ncvar_get(nc,"year")
  }))

  yday <- do.call("c",lapply(filenames,FUN = function(x){
    nc <- ncdf4::nc_open(sprintf("%s/%s",path,x))
    ncdf4::ncvar_get(nc,"day_of_year")
  }))

  # select layers to subset
  layers <- which((years == (year - 1) & yday >= offset) |
                   (years == year & yday < offset))

  # check if all layers are available in the dataset
  # needs a fix for crossing boundaries between decadal subsets !!
  if (length(layers) < 365){
    stop("The selected Berkeley Earth data doesn't cover your time frame!")
  }

  # subset raster data, do this before cropping to make things faster
  delta <- raster::subset(delta, layers)

  # crop data if necessary, no checks yet
  if (!is.null(extent)){
    if(length(extent)==4){
      climatology = raster::crop(climatology, raster::extent(extent))
      delta = raster::crop(delta, raster::extent(extent))
    } else {
      stop("not enough coordinate points to properly crop data!")
    }
  }

  # shift data when offset is < 365
  if (offset < length(layers)){
    doy <- c(offset:length(layers),1:(offset - 1))
  } else {
    doy <- 1:length(layers)
  }

  # calculate absolute temperatures not differences
  # with the mean
  temperature <- raster::stackApply(raster::stack(delta, climatology),
                                   indices = c(doy,1:length(layers)),
                                   fun = sum )
  temperature[temperature == 0] <- NA

  # convert temperature data to matrix
  Ti <- t(raster::as.matrix(temperature))

  # extract georeferencing info to be passed along
  ext <- raster::extent(temperature)
  proj <- raster::projection(temperature)
  size <- dim(temperature)

  cat("calculating daylength \n")

  # grab coordinates
  location <- sp::SpatialPoints(sp::coordinates(temperature),
                               proj4string = sp::CRS(proj))
  location <- t(sp::spTransform(location,
                               sp::CRS("+init=epsg:4326"))@coords[,2:1])

  # create doy vector
  if (offset < length(layers)){
    doy <- c(offset:length(layers),1:(offset - 1))
  } else {
    doy <- 1:length(layers)
  }

  # create daylength matrix
  Li <- lapply(location[1,],
              FUN = function(x){
                daylength(doy = doy, latitude = x)
              })
  Li <- t(do.call("rbind",Li))

  # recreate the validation data structure (new format)
  # but with concatted data
  data <- list("site" = NULL,
              "location" = location,
              "doy" = doy,
              "transition_dates" = NULL,
              "Ti" = Ti,
              "Tmini" = NULL,
              "Tmaxi" = NULL,
              "Li" = Li,
              "Pi" = NULL,
              "VPDi" = NULL,
              "georeferencing" = list("extent" = ext,
                                      "projection" = proj,
                                      "size" = size)
  )

  # assign a class for post-processing
  class(data) <- "phenor_map_data"

  # return the formatted, faster data format
  # either internally or saved as an rds (binary R data file)
  if (internal){
    return(data)
  } else {
    saveRDS(data, file = sprintf("%s/phenor_be_data_%s_%s.rds",path, year))
  }
}

#' Preprocessing of Daymet tiled data into a format which can be ingested
#' by the models in phenor
#'
#' @param path: a path to tiled data
#' @param year: year to process (requires year - 1 to be present)
#' @param tiles: daymet tile number
#' @param offset: offset of the time series in DOY (default = 264, sept 21)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples

# create subset of layers to calculate phenology model output on
format_daymet_tiles = function(path = ".",
                         year = 2014,
                         tiles = 11935,
                         offset = 264,
                         internal = TRUE){

  # some feedback
  cat("calculating average daily temperatures, or load from file \n")

  t1 = sprintf('%s/tmean_%s_%s.tif',path, year - 1, tile)
  t2 = sprintf('%s/tmean_%s_%s.tif',path, year, tile)

  # calculating mean temperature or reading in files
  if (file.exists(t1)) {
    t1 = stack(t1)
  } else {
    t1 = daymet_tmean(
      path = path,
      year = year - 1,
      tile = tile,
      internal = TRUE
    )
  }

  if (file.exists(t2)) {
    t2 = stack(t2)
  } else {
    t2 = daymet_tmean(
      path = path,
      year = year,
      tile = tile,
      internal = TRUE
    )
  }

  cat("subset data \n")
  # create a subset of the data
  t_subset = daymet_subset(stack(t1,t2),
                           offset = offset)
  t_subset_brick = trim(brick(t_subset))

  # convert temperature data to matrix
  Ti = t(as.matrix(t_subset_brick))

  # extract georeferencing info to be passed along
  ext = extent(t_subset_brick)
  proj = projection(t_subset_brick)
  size = dim(t_subset_brick)

  cat("calculating daylength \n")

  # grab coordinates
  location = SpatialPoints(coordinates(t_subset_brick),
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

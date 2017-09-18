#' Preprocessing of Daymet tiled data into a format which can be ingested
#' by the models in phenor
#'
#' @param path: a path to tiled data
#' @param year: year to process (requires year - 1 to be present)
#' @param tiles: daymet tile number
#' @param offset: offset of the time series in DOY (default = 264, sept 21)
#' @return Data file adhering to the phenor modelling input formatting. Data
#' serves as input for model spatial model runs.
#' @keywords phenology, model, preprocessing, climatology
#' @export
#' @examples
#'
#' \dontrun{
#' # run with default settings
#' # looks for daymet average temperature
#' # data in your home directory. These data
#' # are calculated using daymet_tmean() from
#' # the daymetr package
#' daymet_data = format_daymet_tiles()
#'}

# create subset of layers to calculate phenology model output on
format_daymet_tiles = function(path = "~",
                         year = 2014,
                         tile = 11935,
                         offset = 264,
                         internal = TRUE){

  # format paths of the daymet tiles
  tmean_1 = sprintf('%s/tmean_%s_%s.tif',path, year - 1, tile)
  tmean_2 = sprintf('%s/tmean_%s_%s.tif',path, year, tile)
  prcp_1 = sprintf('%s/prcp_%s_%s.nc',path, year - 1, tile)
  prcp_2 = sprintf('%s/prcp_%s_%s.nc',path, year, tile)
  vp_1 = sprintf('%s/vp_%s_%s.nc',path, year - 1, tile)
  vp_2 = sprintf('%s/vp_%s_%s.nc',path, year, tile)

  # process the temperature data
  if (file.exists(tmean_1) & file.exists(tmean_2) ) {
    tmean_1 = stack(tmean_1); tmean_2 = stack(tmean_2)
    tmean_subset = daymet_subset(stack(tmean_1,tmean_2), offset = offset)
    tmean_subset_brick = trim(brick(tmean_subset))
    Ti = t(raster::as.matrix(tmean_subset_brick))
  } else {
    stop("No average daily temperature files are found\n
         please generate these files first using daymet_tmean()\n
         from the daymetr package.")
  }

  # process the precipitation data
  if(file.exists(p1) & file.exists(p2)){
    prcp_1 = stack(prcp_1); prcp_2 = stack(prcp_2)
    prcp_subset = daymet_subset(stack(prcp_1,prcp_2), offset = offset)
    prcp_subset_brick = trim(brick(prcp_subset))
    Pi = t(raster::as.matrix(prcp_subset_brick))
  } else {
    warning("Correct precipitation files are not provided, will return NULL.")
    Pi = NULL
  }

  # process the VPD data
  if(file.exists(vp_1) & file.exists(vp_2)){
    vp_1 = stack(vp_1); vp_2 = stack(vp_2)
    vp_subset = daymet_subset(stack(vp_1,vp_2), offset = offset)
    vp_subset_brick = trim(brick(vp_subset))
    VPDi = t(raster::as.matrix(vp_subset_brick))
  } else {
    warning("Correct vapour pressure files are not provided, will return NULL.")
    VPDi = NULL
  }

  # extract georeferencing info to be passed along
  ext = extent(t_subset_brick)
  proj = projection(t_subset_brick)
  size = dim(t_subset_brick)

  # grab coordinates
  location = SpatialPoints(coordinates(t_subset_brick),
                           proj4string = CRS(proj))
  location = t(spTransform(location, CRS("+init=epsg:4326"))@coords[,2:1])

  # create doy vector
  if (offset < 365){
    doy = c(offset:365,1:(offset - 1))
  } else {
    doy = 1:365
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
              "Pi" = Pi,
              "VPDi" = VPDi,
              "altitude" = NULL,
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

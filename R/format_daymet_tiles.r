#' Preprocessing of Daymet tiled data into a format which can be ingested
#' by the models in phenor
#'
#' @param path a path to tiled data
#' @param year year to process (requires year - 1 to be present)
#' @param tile daymet tile number
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param internal return results as an R variable or save as an RDS file
#' @return Data file adhering to the phenor modelling input formatting. Data
#' serves as input for model spatial model runs. The only required data input
#' is the mean daily temperature, precipitation and VPD can be included but
#' are not mandatory. Depending on the models used, including only temperature
#' will keep file sizes down and processing time lower.
#' @keywords phenology, model, preprocessing, climatology
#' @export
#' @examples
#'
#' # run with default settings
#' # looks for daymet average temperature
#' # data in your home directory. These data
#' # are calculated using daymet_tmean() from
#' # the daymetr package
#' \dontrun{
#' daymet_data = format_daymet_tiles()
#'}

# create subset of layers to calculate phenology model output on
format_daymet_tiles = function(path = tempdir(),
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

  # MOVE INTO A FUNCTION !!
  # process the temperature data
  if (file.exists(tmean_1) & file.exists(tmean_2) ) {
    tmean_1 = raster::stack(tmean_1)
    tmean_2 = raster::stack(tmean_2)
    tmean_subset = daymetr::daymet_grid_offset(raster::stack(tmean_1,tmean_2),
                                 offset = offset)
    tmean_subset_brick = raster::trim(raster::brick(tmean_subset))
    Ti = t(raster::as.matrix(tmean_subset_brick))
  } else {
    stop("No average daily temperature files are found\n
         please generate these files first using daymet_tmean()\n
         from the daymetr package.")
  }

  # process the precipitation data
  if(file.exists(prcp_1) & file.exists(prcp_2)){
    prcp_1 = raster::stack(prcp_1)
    prcp_2 = raster::stack(prcp_2)
    prcp_subset = daymetr::daymet_grid_offset(raster::stack(prcp_1,prcp_2),
                                offset = offset)
    prcp_subset_brick = raster::trim(raster::brick(prcp_subset))
    Pi = t(raster::as.matrix(prcp_subset_brick))
  } else {
    warning("Correct precipitation files are not provided, will return NULL.")
    Pi = NULL
  }

  # process the VPD data
  if(file.exists(vp_1) & file.exists(vp_2)){
    vp_1 = raster::stack(vp_1)
    vp_2 = raster::stack(vp_2)
    vp_subset = daymetr::daymet_grid_offset(raster::stack(vp_1,vp_2),
                              offset = offset)
    vp_subset_brick = raster::trim(raster::brick(vp_subset))
    VPDi = t(raster::as.matrix(vp_subset_brick))
  } else {
    warning("Correct vapour pressure files are not provided, will return NULL.")
    VPDi = NULL
  }

  # extract georeferencing info to be passed along
  ext = raster::extent(tmean_subset_brick)
  proj = raster::projection(tmean_subset_brick)
  size = dim(tmean_subset_brick)

  # grab coordinates
  location = sp::SpatialPoints(sp::coordinates(tmean_subset_brick),
                           proj4string = sp::CRS(proj))
  location = t(sp::spTransform(location,
                               sp::CRS("+init=epsg:4326"))@coords[,2:1])

  # create doy vector
  if (offset < 365){
    doy = c(offset:365,1:(offset - 1))
  } else {
    doy = 1:365
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
              "VPDi" = VPDi,
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
    saveRDS(data, file = sprintf("%s/phenor_daymet_data_%s_%s.rds",path, year, tile))
  }
}

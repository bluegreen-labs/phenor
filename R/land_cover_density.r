#' Calculates land cover density values based upon a MODIS MCD12Q1
#' land cover input map (src_raster) and a coarser destination raster.
#'
#' The algorithm counts the number of occurences in a given pixel of
#' the destination raster. Both maps should be in EPSG:4326.
#'
#' @param dest_raster area to calculate the statistics for, with a coarser
#' resolution than the 500m MCD12Q1 data
#' @param src_raster a MCD12Q1 map in epsg:4326, if NULL the internal
#' CONUS map will be used. This maps is produced by calculating the most
#' common land cover class between 2001 - 2009 for CONUS.
#' @param lc_classes land cover classes to calculate the statistics for.
#' only IGBP classes (1 - 16) are supported. Takes a vector e.g. c(1, 4)
#' @param path path to output your generated data to if not outputting
#' it to the R console / workspace, default = "~"
#' @param internal TRUE / FALSE, if true no files are written to disk and
#' a raster stack is returned to the R command line
#' @keywords phenology, model, validation, comparison
#' @export
#' @examples
#'
#' # will return a land cover density map for
#' # the evergreen needleleaf class (1) and
#' # the deciduous broadleaf classs (4) in layers
#' # 1 and 2 of the raster stack lc_dens
#' \dontrun{
#' lc_dens = land_cover_density(lc_classes = c(1, 4))
#' }

land_cover_density = function(src_raster = NULL,
                              dest_raster = NULL,
                              lc_classes = c(1,4,5,10),
                              path = "~",
                              internal = FALSE){

  # MCD12Q1 land cover class file
  if(is.null(src_raster)){
    stop("No source raster provided.")
  } else {
    if(attr(class(src_raster), "package") == "raster"){
      lc = src_raster
    } else {
      lc = raster::raster(src_raster)
    }
  }

  # destination raster
  if(is.null(dest_raster)){
    stop("No destination raster provided.")
  } else {
    if(attr(class(src_raster), "package") == "raster"){
      dest_raster = src_raster
    } else {
      dest_raster = raster::raster(dest_raster)
    }
  }

  # set projection (lat long)
  lat_lon = sp::CRS("+init=epsg:4326")

  # extract size info from the destination raster
  rows = nrow(dest_raster)
  cols = ncol(dest_raster)

  # number the pixels for further reference
  zones = dest_raster
  zones[] = 1:prod(rows,cols)

  # now loop over all land cover classes you want
  # to extract from the original map for the grid
  # cells you set by the dest_raster
  for (i in lc_classes){

    # get cell numbers for all pixels of class i
    cell_v = raster::Which(lc == i, cells = TRUE)

    # extract coordinates for all pixels of class i
    xy = raster::xyFromCell(lc, cell = cell_v)

    # read in the coordinates and assign them a projection
    ll = sp::SpatialPoints(xy, lat_lon)

    # get the zone number for pixels with location ll
    zone_numbers = raster::extract(zones, ll)

    # merge results, subset based upon pixel class i
    subs_matrix = data.frame(zone_numbers,
                             rep(1, length(cell_v)))

    # count the number of forest pixels in a given cell
    # of the dest_raster format
    count_v = by(data = subs_matrix[,2],
                 INDICES = subs_matrix[,1],
                 FUN = function(x,...){length(x)},
                 na.rm = TRUE)

    # get the maximum count
    max_v = max(count_v)

    # get cell ids
    id = as.numeric(names(count_v))

    # define substitution matrix
    smatrix = data.frame(as.vector(id),
                         as.vector(count_v))

    # substitute values
    tmp_cover = raster::subs(zones, smatrix, which = 2)

    # normalize
    tmp_cover = tmp_cover / max_v

    if (i != lc_classes[1]){
      lc_cover = tmp_coverage
    } else {
      lc_cover = raster::stack(lc_cover, tmp_cover)
    }
  }

  # assign names to the layers
  names(lc_cover) = lc_classes

  if(internal){
    # return the data as a raster stack
    # beware when using the internal option
    # this file might get really big depending
    # on the level op upscaling
    return(lc_cover)
  } else {
    # write everything to file
    raster::writeRaster(lc_cover,
                        sprintf('%s/igbp.tif',path),
                        overwrite = TRUE,
                        options = c("COMPRESS=DEFLATE"))
  }
}

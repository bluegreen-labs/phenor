#' Calculates land cover density
#'
#' Use MODIS MCD12Q1 land cover input map (src_raster) and a coarser
#' destination raster (E-OBS, CMIP5) to calculate land cover representation
#' of these coarser cells.
#'
#' The algorithm counts the number of occurences in a given pixel of
#' the destination raster. Both maps should be in EPSG:4326.
#'
#' @param lc_raster a MCD12Q1 map in epsg:4326.
#' @param dest_raster area to calculate the statistics for, with a coarser
#' resolution than the 500m MCD12Q1 data
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
#' \dontrun{
#' # will return a land cover density map for
#' # the evergreen needleleaf class (1) and
#' # the deciduous broadleaf classs (4)
#' lc_dens = pr_calc_land_cover(lc_classes = c(1, 4),
#'                              lc_raster = "~/MCD12Q1_igbp_map.tif"
#'                              dest_raster = "~/1_degree_lat_lon_map.tif")
#' }

pr_calc_land_cover = function(
  lc_raster,
  dest_raster,
  lc_classes = c(1,4,5,10),
  path = tempdir(),
  internal = FALSE
){

  # MCD12Q1 land cover class file
  if(missing(lc_raster)){
    stop("No source raster provided.")
  } else {
    if(!class(lc_raster) != "RasterLayer"){
      if(is.character(lc_raster)){
        if(!file.exists(deparse(substitute(lc_raster)))){
          stop("No source raster provided.")
        }
      }
    } else {
        lc_raster = raster::raster(lc_raster)
    }
  }

  # destination raster
  if(missing(dest_raster)){
    stop("No source raster provided.")
  } else {
    if(!class(dest_raster) != "RasterLayer"){
      if(is.character(dest_raster)){
        if(!file.exists(deparse(substitute(dest_raster)))){
          stop("No source raster provided.")
        }
      }
    } else {
      dest_raster = raster::raster(dest_raster)
    }
  }

  # clean up existing files, should some remain in tempdir
  # due to interupted processing
  img_files <- list.files(tempdir(),
                          utils::glob2rx("fractional_land_cover_*.tif"),
                          full.names = TRUE)
  if(length(img_files)!=0) {
    file.remove(img_files)
  }

  # extract size info from the destination raster
  rows = nrow(dest_raster)
  cols = ncol(dest_raster)

  # number the pixels for further reference
  zones = dest_raster
  zones[] = 1:prod(rows,cols)

  # now loop over all land cover classes you want
  # to extract from the original map for the grid
  # cells you set by the dest_raster
  lapply(lc_classes, function(i){

    # feedback
    cat(sprintf("processing land cover class: %s \n", i))

    # get cell numbers for all pixels of class i
    cell_v = raster::Which(lc_raster == i, cells = TRUE)

    # extract coordinates for all pixels of class i
    xy = raster::xyFromCell(lc_raster, cell = cell_v)

    # read in the coordinates and assign them a projection
    ll = sp::SpatialPoints(xy, sp::CRS("+init=epsg:4326"))

    # get the zone number for pixels with location ll
    zone_numbers = raster::extract(zones, ll)

    # merge results, subset based upon pixel class i
    subs_matrix = data.frame(zone_numbers,
                             rep(1, length(cell_v)))

    # count the number of xyz pixels in a given cell
    # of the dest_raster format
    count_v = by(data = subs_matrix[,2],
                 INDICES = subs_matrix[,1],
                 FUN = function(x,...){length(x)},
                 na.rm = TRUE)

    # get the maximum count
    max_v = max(count_v, na.rm = TRUE)

    # get cell ids
    id = as.numeric(names(count_v))

    # define substitution matrix
    smatrix = data.frame(as.vector(id),
                         as.vector(count_v))

    # substitute values
    tmp_cover = raster::subs(zones, smatrix, which = 2)

    # normalize
    tmp_cover = tmp_cover/max_v

    # store data in a temporary tif file
    # (doesn't take up memory which needs to be cleared for efficiency)
    # output data
    raster::writeRaster(tmp_cover,
                sprintf("%s/fractional_land_cover_%02d.tif", tempdir(), i),
                overwrite = TRUE)

    # clear junk from memory
    gc()
  })

  # list tif files
  fractional_cover = raster::stack(
    list.files(tempdir(),
               utils::glob2rx("fractional_land_cover_*.tif"),
               full.names = TRUE))

  # assign names to the layers
  names(fractional_cover) = paste0("igbp_", lc_classes)

  if(internal){
    # return the data as a raster stack
    return(fractional_cover)
  } else {
    # write everything to file
    raster::writeRaster(fractional_cover,
                        sprintf('%s/fractional_land_cover.tif',path),
                        overwrite = TRUE,
                        options = c("COMPRESS=DEFLATE"))
  }
}

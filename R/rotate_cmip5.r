#' Rotate CMIP5 NetCDF data cubes
#'
#' Rotate NetCDF files faster than the default
#' raster rotate() command, as data is cropped
#' before any translations (reducing memory load)

#' @param r raster layer, stack or brick
#' @param extent vector with coordinates defining the region of interest defined
#' as xmin, xmax, ymin, ymax in lat/lon (default = c(-74,-65, 40, 48))
#' @return Returns raster data in the same format (layer, stack, brick) as
#' the input, but rotated from a climatologica longitude (0 - 360) format to
#' a normal longitude format (-180 to 180).
#' @export

rotate_cmip5 <- function(r = NULL,
                           extent = c(-74,-65, 40, 48)){

  # duplicate the extent values
  # to swap out the longitude ones
  extent2 = extent

  # check minimum longtiude value
  if (extent[1] < 0){
    extent2[1] = 360 + extent[1]
  } else {
    m = crop(r, extent)
  }

  # check maximum longitude value
  if (extent[2] < 0){
    extent2[2] = 360 + extent[2]
    m = shift(crop(r, extent2), x = -360)
  } else {
    extent2[2] = 360
    extent[1] = 0
    m = mosaic(shift(crop(r, extent2), x = -360),
               crop(r, extent),fun=max)
  }

  # return mosaic
  return(m)
}

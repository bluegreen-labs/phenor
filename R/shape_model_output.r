#' (re)shape model output based upon the class of the input data
#' and valid model estimates. Mainly, reshapes data to a spatial
#' raster format when required.
#'
#' @param data input data generated using the format_*() functions
#' @param doy phenophase estimates as a doy value
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' # Should not be run stand-alone
#' \dontrun{
#' shape_model_output(data = data, doy = doy)
#'}

shape_model_output = function(data, doy){
  if(class(data) == "phenor_map_data"){
    r = raster::raster(nrows = data$georeferencing$size[1],
                       ncols = data$georeferencing$size[2])
    raster::extent(r) = data$georeferencing$extent
    sp::proj4string(r) = sp::CRS(data$georeferencing$projection)
    r[] = as.vector(doy)
    r[r==9999] = NA
    return(r)
  } else {
    return(doy)
  }
}

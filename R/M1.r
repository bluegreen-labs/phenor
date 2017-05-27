#' M1 model as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data: a nested list of data with on location:
#' 1. the date (doy or long format)
#' 2. the temperature data
#' 3. the photoperiod data (NA when not needed)
#' 4. a vector or matrix with necessary constants (NA when not needed)
#'    - long term mean temperature
#'    - latitude
#'    - etc...
#' @param par: a vector of parameter values, this is functions specific
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = TT(data = data, par = par)
#'}

M1 = function(par, data, plot = FALSE){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 = par[1]
  T_base = par[2]
  k = par[3]
  F_crit = par[4]

  # create forcing/chilling rate vector
  # forcing
  Rf = data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf = ((data$Li / 24) ^ k) * Rf
  Rf[1:t0,] = 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = data$doy[which(cumsum(xt) >= F_crit)[1]]
    doy[is.na(doy)] = 9999
    return(doy)
  })

  # set export format, either a rasterLayer
  # or a vector
  if(is.null(data$site)){
    r = raster(nrows = data$georeferencing$size[1],
               ncols = data$georeferencing$size[2])
    extent(r) = data$georeferencing$extent
    proj4string(r) = CRS(data$georeferencing$projection)
    r[] = doy
    r[r==9999] = NA
    return(r)
  } else {
    return(doy)
  }
}

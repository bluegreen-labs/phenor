#' Thermal Time model as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data: a nested list of data
#' @param par: a vector of parameter values, this is functions specific
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = TT(data = data, par = par)
#'}

TT = function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 = round(par[1])
  T_base = par[2]
  F_crit = par[3]

  # simple degree day sum setup
  Rf = data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf[1:t0,] = 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = data$doy[which(cumsum(xt) >= F_crit)[1]]
    doy[is.na(doy)] = 9999
    return(doy)
  })

  # set export
  if(is.null(data$site)){
    r = raster(nrows = data$size[1], ncols = data$size[2])
    extent(r) = data$extent
    proj4string(r) = CRS(data$projection)
    r[] = doy
    r[r==9999] = NA
    return(r)
  } else {
    return(doy)
  }
}

#' Parallel model as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' using a bell shaped chilling response
#'
#' @param data: a nested list of data
#' @param par: a vector of parameter values, this is functions specific
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = PAb(data = data, par = par)
#'}

PAb = function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 9){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 = round(par[1])
  t0_chill = round(par[2])
  T_base = par[3]
  C_a = par[4]
  C_b = par[5]
  C_c = par[6]
  C_ini = par[7]
  F_crit = par[8]
  C_req = par[9]

  # chilling
  Rc = 1 / ( 1 + exp( C_a * (data$Ti - C_c) ^ 2 + C_b * (data$Ti - C_c) ))
  Rc[1:t0_chill,] = 0
  Sc = apply(Rc,2, cumsum)

  # chilling requirement has to be met before
  # accumulation starts (binary choice)
  k = C_ini + Sc * (1 - C_ini)/C_req
  k[Sc >= C_req] = 1

  # forcing
  Rf = data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf = Rf * k
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

#' Accumulated growing season index model
#' as defined by Xin et al. 2015 (Rem. Sens. Env.)
#'
#' @param data: input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par: a vector of parameter values, this is functions specific
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = AGSI(data = data, par = par)
#'}

AGSI = function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 2){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 = round(par[1])
  F_crit = par[2]

  # set constants as defined by
  # jolly et al. 2005 (Glob. Change. Biol.)
  # these are physiological limits
  Tmmin = -2
  Tmmax = 5
  VPDmin = 900
  VPDmax = 4100
  photo_min = 10
  photo_max = 11

  # calculate matrices
  Tmin = (data$Tmini - Tmmin)/(Tmmax - Tmmin)
  VPD = (data$VPDi - VPDmin)/(VPDmax - VPDmin)
  photo = (data$Li - photo_min)/(photo_max - photo_min)

  # criteria
  Tmin[which(Tmin <= Tmmin)] = 0
  Tmin[which(Tmin >= Tmmax)] = 1

  VPD[which(VPD >= VPDmax)] = 0
  VPD[which(VPD <= VPDmin)] = 1

  photo[which(photo <= photo_min)] = 0
  photo[which(photo >= photo_max)] = 0

  # calculate the index for every day
  GSI = Tmin * VPD * photo
  GSI[1:t0,] = 0

  # DOY of budburst criterium
  doy = apply(GSI,2, function(xt){
    doy = data$doy[which(cumsum(xt) >= F_crit)[1]]
    doy[is.na(doy)] = 9999
    return(doy)
  })

  # set export format, either a rasterLayer
  # or a vector
  if(class(data) == "phenor_map_data"){
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

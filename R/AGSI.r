#' Accumulated growing season index model (1)
#' as defined by Xin et al. 2015 (Rem. Sens. Env.)
#'
#' The starting point of accumulation is not clearly indicated
#' in the publication so we assume December 21.
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = AGSI(data = data, par = par)
#'}

AGSI = function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop(sprintf("model parameter(s) out of range (too many, too few, %s provided)",
                 length(par)))
  }

  # exit the routine if data is missing
  if (is.null(data$Tmini) |
      is.null(data$Tmaxi) |
      is.null(data$VPDi)  |
      is.null(data$Li)
      ){
    stop("Not all required driver data is available")
  }

  # set start of accumulation period
  t0 = which(data$doy == -11)
  if (length(t0)==0){
    t0 = 1
  }

  # extract the parameter values from the
  # par argument for readability
  Tmmin = par[1]
  Tmmax = par[2]
  F_crit = par[3]
  VPDmin = par[4]
  VPDmax = par[5]
  photo_min = par[6]
  photo_max = par[7]

  # rescaling all parameters between 0 - 1 for the
  # acceptable ranges, if outside these ranges
  # set to 1 (unity) or 0 respectively
  Tmin = (data$Tmini - Tmmin)/(Tmmax - Tmmin)
  VPD = (data$VPDi - VPDmin)/(VPDmax - VPDmin)
  photo = (data$Li - photo_min)/(photo_max - photo_min)

  # set outliers to 1 or 0
  Tmin[which(data$Tmini <= Tmmin)] = 0
  Tmin[which(data$Tmini >= Tmmax)] = 1

  VPD[which(data$VPDi >= VPDmax)] = 0
  VPD[which(data$VPDi <= VPDmin)] = 1

  photo[which(data$Li <= photo_min)] = 0
  photo[which(data$Li >= photo_max)] = 1

  # calculate the index for every value
  # but set values for time t0 to 0
  GSI = Tmin * VPD * photo
  GSI[1:t0,] = 0

  # DOY of budburst criterium as calculated
  # by cummulating the GSI hence AGSI
  doy = apply(GSI,2, function(xt){
    doy = data$doy[which(cumsum(xt) >= F_crit)[1]]
    doy[is.na(doy)] = 9999
    return(doy)
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

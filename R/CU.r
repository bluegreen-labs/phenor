#' CU chilling degree model
#' Chen et al. 2017 (Journal of Ecology)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model, tropical forest
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = CU(data = data, par = par)
#'}

CU <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  ni = round(par[1])
  nd = round(par[2])
  T_base = par[3]
  F_crit = par[4]

  # simple chilling degree day sum setup
  # with lagged response
  data$Ti[data$Ti > T_base] = 0
  data$Ti <- abs(data$Ti)

  # calculate floral induction time series
  Fd <- apply(data$Ti, 2, function(x){
    Fi <- rollapply(x, ni, sum, align = "right", fill = 0, na.rm = TRUE)
    Fd <- c(rep(0,nd),Fi[1:(length(Fi) - nd)])
  })

  # DOY of budburst criterium
  doy = apply(Fd, 2, function(xt){
    doy = data$doy[which(xt >= F_crit)[1]]
    doy[is.na(doy)] = 9999
    return(doy)
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' DormPhot model as defined in
#' Caffarra, Donnelly and Chuine 2011 (Clim. Res.)
#' parameter ranges are taken from Basler et al. 2016
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
#' estimate = DP(data = data, par = par)
#'}

DP <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 11){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  a = par[1]
  b = par[2]
  c = par[3]
  d = par[4]
  e = par[5]
  f = par[6]
  g = par[7]
  F_crit = par[8]
  C_crit = par[9]
  L_crit = par[10]
  D_crit = par[11]

  # set the t0 value if necessary (~Sept. 1)
  t0 = which(data$doy < 1 & data$doy == -121)

  # dormancy induction
  # (this is all vectorized doing cell by cell multiplications with
  # the sub matrices in the nested list)
  DR = 1/(1 + exp(a * (data$Ti - b))) * 1/(1 + exp(10 * (data$Li - L_crit)))
  if (!length(t0) == 0){
    DR[1:t0,] = 0
  }
  DS = apply(DR,2, cumsum)

  # chilling
  CR = 1/(1 + exp(c * (data$Ti - d)^2 + (data$Ti - d) ))
  CR[DS < D_crit] = 0
  CS = apply(CR,2, cumsum)

  # forcing
  dl50 = 24 / (1 + exp(f * (CS - C_crit))) # f >= 0
  t50 = 60 / (1 + exp(g * (data$Li - dl50))) # g >= 0
  Rf = 1/(1 + exp(e * (data$Ti - t50))) # e <= 0
  Rf[DS < D_crit] = 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = data$doy[which(cumsum(xt) >= F_crit)[1]]
    doy[is.na(doy)] = 9999
    return(doy)
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

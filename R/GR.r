#' Thermal Time grassland model as defined in
#' Garcia-Mozo et al. 2009 (Agr. For. Meteorlogy)
#'
#' @param data: input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par: a vector of parameter values, this is functions specific
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = GR(data = data, par = par)
#'}

GR = function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  b = par[1]
  c = par[2]
  F_crit = par[3]
  P_crit = par[4]
  L_crit = par[5]

  # find day on which in the past 7 days more than
  # P_crit rain fell
  P_star = ifelse(data$Li >= L_crit, 1, 0)

  # rainfall accumulation
  R_star = data$Pi * P_star
  R_star = rollapply(R_star,
                width = 7,
                FUN = sum,
                align = "left",
                fill = 0,
                by.column = TRUE)

  # find day on which in the past 7 days more than
  # P_crit rain fell
  R_star = apply(R_star, 2, function(xt) {
    v = rep(0, length(xt))
    l = which(xt >= P_crit)[1]
    if (!is.na(l)) {
      v[l:length(xt)] = 1
    }
    return(v)
  })

  # create forcing/chilling rate vector
  # forcing
  Rf = 1 / (1 + exp(-b * (data$Ti - c))) * R_star

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = data$doy[which(cumsum(xt) >= F_crit)[1]]
    doy[is.na(doy)] = 9999
    return(doy)
  })
  return(doy)
}

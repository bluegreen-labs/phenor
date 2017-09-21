#' Thermal Time grassland model which includes a pulse response
#' precipitatoin trigger as defined in
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

  # THIS FUNCTION IS TOO SLOW, ONLY CALCULATE UNTIL
  # THRESHOLD IS REACHED.
  data$Pi = data$Pi * P_star
  R_star = matrix(0,nrow(data$Pi),ncol(data$Pi))

  # This avoids inefficient moving window
  # approaches which are computationally
  # expensive (break the routine when the
  # criterium is met instead of a full run
  # for a year)
  for (i in 1:ncol(R_star)){
    for (j in 1:nrow(R_star)){
      if (j == nrow(R_star) - 7){
        break
      }
      if(sum(data$Pi[j:(j+7),i]) >= P_crit){
        R_star[j:nrow(R_star),i] = 1
        break
      }
    }
  }

  # create forcing/chilling rate vector
  # forcing
  Rf = 1 / (1 + exp(b * (data$Ti + c))) * R_star

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
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

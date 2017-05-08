#' Linear model as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data: a nested list of data
#' @param par: a vector of parameter values, this is functions specific
#' @param spring: a vector defining spring as a couple of DOY values
#' default is March, April, May or DOY 60 - 151
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = LIN(data = data, par = par)
#'}

LIN = function(par, data, spring = c(60,151)){

  # exit the routine as some parameters are missing
  if (length(par) != 2){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  a = par[1]
  b = par[2]

  # calculate location of "spring" values
  spring_loc = data$doy %in% seq(spring[1],spring[2])

  # calculate the mean temperature for this range
  mean_spring_temp = apply(data$Ti,2,function(xt){
    mean(xt[spring_loc])
  })

  # return fit values
  return( a * mean_spring_temp + b )
}

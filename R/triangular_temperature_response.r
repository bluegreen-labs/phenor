#' Triangular temperature response function as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param T a vector or matrix of temperatures
#' @param T_opt optimal temperature
#' @param T_min minimum viable temperature
#' @param T_max maximum viable temperature
#' @keywords phenology, model, temperature response
#' @export
#' @examples
#'
#' T_response = triangular_temperature_response(T = 0:45)
#' \dontrun{
#' plot(0:45, T_response, type = "l")
#'}

triangular_temperature_response <- function(T = -10:45,
                                    T_opt = 10,
                                    T_min = 1,
                                    T_max = 15){

  # sanity checks
  if (T_opt >= T_max || T_opt <= T_min || T_max <= T_min){
    T[] = NA
    return(T)
  }

  # find locations of rising and falling
  # part of the triangular function
  loc_rising <- which(T < T_opt & T >= T_min)
  loc_falling <- which(T <= T_max & T >= T_opt)

  # fill this vector according to a triangular
  # ruleset function

  # set out of range values to 0
  T[T < T_min | T > T_max] = 0

  # convert temperatures
  T[loc_rising] <- (T[loc_rising] - T_min)/(T_opt - T_min)
  T[loc_falling] <- 1 - (T[loc_falling] - T_opt)/(T_max - T_opt)

  # returns the converted temperature data
  return(T)
}

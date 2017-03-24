#' Calculates day length (in hours) and the solar elevation above the
#' ecliptic plane based upon latitude and a day of year value.
#' 
#' @param doy: a vector with doy values 1 - 365(6)
#' @param latitude: a given latitude
#' @keywords solar, ephemerids
#' @export
#' @examples
#' ephem <- daylength(1,51)
#' 
#' # hours of sunlight on day 1 of the year
#' print(ephem)

daylength = function(doy,latitude){
  
  # define constant
  p = 0
  
  # Forsythe et al. 1995 eq. 1
  Omega = 0.2163108 + 2 * atan(0.9671396 * tan (0.00860 * (doy - 186)))
  
  # eq. 2
  Phi = asin(0.39795* cos(Omega))
  
  # eq. 3 / returns daylength D
  DL = 24 - 24 / pi * acos((sin(p*pi/180)+sin(latitude*pi/180)*sin(Phi))/(cos(latitude*pi/180)*cos(Phi)))
  
  for (i in 1:length(DL)){
    l <- DL[i-1] > 20
    if (length(l)==0){
      l <- FALSE
    }
    l[is.na(l)] <- FALSE
    if ( l  & is.na(DL[i]) ){
      DL[i] <- 24
    }
    if (is.na(DL[i])){
      DL[i-1] <- 0
    }
  }
  return(list(DL,Phi*180/pi))
}
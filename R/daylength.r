#' Calculates day length (in hours) and the solar elevation above the
#' ecliptic plane based upon latitude and a day of year value.
#' This routine uses Forsythe et al. 1995.
#'
#' @param doy: a vector with doy values 1 - 365(6)
#' @param latitude: a given latitude
#' @keywords solar, ephemerids
#' @export
#' @examples
#'
#' \dontrun{
#' # calcualte the hours of sunlight and solar elevation on day of year 1
#' ephem = daylength(1,51)
#' print(ephem)
#' }

daylength = function(doy, latitude) {

  # convert to numeric to be sure
  latitude = as.numeric(latitude)

  # define constant
  p = 0

  # degrees to radial conversion
  conv = pi / 180

  # Forsythe et al. 1995 eq. 1
  Omega = 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))

  # eq. 2
  Phi = asin(0.39795 * cos(Omega))

  # eq. 3 / returns daylength D
  DL = suppressWarnings(24 - 24 / pi * acos((sin(p * conv) + sin(latitude * conv) *
                              sin(Phi)) / (cos(latitude * conv) * cos(Phi))))

  # convert declination to solar elevation (above ecliptica) or
  # 90 - zenith angle
  solar_elev = 90 - acos(sin(latitude * conv) * sin(Phi) + cos(latitude * conv) * cos(Phi)) * 180 / pi

  # addresses NA / inf values in output due to the asymptotic nature of the
  # daylength function
  if (latitude > 0){
  for (i in 1:length(DL)) {

    l <- DL[i - 1] > 10

    if (length(l) == 0) {
      l <- FALSE
    }

    l[is.na(l)] <- FALSE

    if (l  & is.na(DL[i])) {
      DL[i] <- 24
    }
    if (is.na(DL[i])) {
      DL[i - 1] <- 0
    }
  }
  } else {
    for (i in 1:length(DL)) {

      l <- DL[i - 1] < 10

      if (length(l) == 0) {
        l <- FALSE
      }

      l[is.na(l)] <- FALSE

      if (l  & is.na(DL[i])) {
        DL[i] <- 0
      }
      if (is.na(DL[i])) {
        DL[i - 1] <- 24
      }
    }
  }
  return(list(DL, solar_elev))
}

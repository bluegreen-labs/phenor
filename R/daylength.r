#' Calculates day length (in hours) and the solar elevation above the
#' ecliptic plane based upon latitude and a day of year, and year values.
#' according to H.Glarner (http://herbert.gandraxa.com/length_of_day.xml)
#'
#' Due to heterogeneous date formats with years normalized to 365 days
#' I do not apply a leap year correction.
#'
#' @param doy: a vector with doy values 1 - 365(6)
#' @param latitude: a given latitude
#' @keywords solar, ephemerids
#' @export
#' @examples
#'
#' \dontrun{
#' # calcualte the hours of sunlight and solar elevation on day of year 1
#' length_of_day = daylength(1, 51, 2000)
#' print(length_of_day)
#' }

daylength = function(doy, latitude) {

  # set constants
  latitude = (pi / 180) * latitude

  # Correct for winter solistice
  doy = doy + 11

  # earths ecliptic
  j = pi / 182.625
  axis = (pi / 180) * 23.439

  # calculate daylength for all days
  dl = lapply(doy, function(x){

    #Exposed radius part between sun's zenith and sun's circle
    m = 1 - tan(latitude) * tan(axis * cos(j * x))

    # sun never appears or disappears
    if (m < 0) { m = 0 }
    if (m > 2) { m = 2 }

    # Exposed fraction of the sun's circle
    b = acos(1 - m) / pi

    # Daylength (lat,day)
    return(b * 24)
  })

  # return the daylength vector
  return(unlist(dl))
}

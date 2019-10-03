#' Calculates AICc values for a set of measured and predicted values
#' together with the number of model parameters used
#'
#' @param measured a vector with measurement values to smooth
#' @param predicted a vector with dates / time steps
#' @param k optional values to weigh the loess fit with
#' @return returns the AIC for measured and predicted values
#' for use in model comparison and selection
#' @keywords model selection, AIC, Akaike's Information Criterion
#' @export
#' @examples
#'
#' \dontrun{
#'
#' model_AIC <- AICc(measured, predicted, k)
#'
#' }

# custom AIC function which accepts loess regressions
AICc <- function(measured, predicted, k){

  # calculate number of observations
  n <- length(measured)

  # calculatue residual sum of squares
  RSS <- sum((measured - predicted)^2)

  # AIC
  AIC <- 2*k + n * log(RSS/n)

  # AICc
  AICc <- AIC + (2 * k * (k + 1)) / (n - k - 1)

  # return both AIC
  return(list("AIC" = AIC,
              "AICc" = AICc))
}

#' Rotate CMIP5 NetCDF data cubes
#'
#' Rotate NetCDF files faster than the default
#' raster rotate() command, as data is cropped
#' before any translations (reducing memory load)
#'
#' @param r CMIP5 raster layer, stack or brick
#' @param extent vector with coordinates defining the region of interest defined
#' as xmin, xmax, ymin, ymax in lat/lon (default = c(-74,-65, 40, 48))
#' @return Returns raster data in the same format (layer, stack, brick) as
#' the input, but rotated from a climatologica longitude (0 - 360) format to
#' a normal longitude format (-180 to 180).
#' @export

rotate_cmip5 <- function(r = NULL,
                         extent = c(-74, -65, 40, 48)){

  # duplicate the extent values
  # to swap out the longitude ones
  extent2 = extent

  # check minimum longtiude value
  # correct for the offset
  if (extent[1] < 0){
    extent2[1] = 360 + extent[1]
  }

  # check maximum longitude value
  # check for offset
  if (extent[2] < 0){
    extent2[2] = 360 + extent[2]
    m = raster::shift(raster::crop(r, extent2), -360)
  } else if (extent[1] < 0 & extent[2] > 0) {
    extent2[2] = 360
    extent[1] = 0
    m = raster::mosaic(raster::shift(raster::crop(r, extent2), -360),
                       raster::crop(r, extent),fun=max)
  } else {
    m = raster::crop(r, extent)
  }

  # return mosaic
  return(m)
}

#' (re)shape model output based upon the class of the input data
#' and valid model estimates. Mainly, reshapes data to a spatial
#' raster format when required.
#'
#' @param data input data generated using the format_*() functions
#' @param doy phenophase estimates as a doy value
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' # Should not be run stand-alone
#' \dontrun{
#' shape_model_output(data = data, doy = doy)
#'}

shape_model_output <- function(data, doy){
  if(class(data) == "phenor_map_data"){
    r = raster::raster(nrows = data$georeferencing$size[1],
                       ncols = data$georeferencing$size[2])
    raster::extent(r) <- data$georeferencing$extent
    sp::proj4string(r) <- sp::CRS(data$georeferencing$projection)
    r[] <- as.vector(doy)
    r[ r == 9999 ] <- NA
    return(r)
  } else {
    return(doy)
  }
}

#' Calculates day length (in hours) and the solar elevation above the
#' ecliptic plane based upon latitude and a day of year, and year values.
#' according to H.Glarner (http://herbert.gandraxa.com/length_of_day.xml)
#'
#' Due to heterogeneous date formats with years normalized to 365 days
#' I do not apply a leap year correction.
#'
#' @param doy a vector with doy values 1 - 365(6)
#' @param latitude a given latitude
#' @return a daylength vector
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
  latitude <- (pi / 180) * latitude

  # Correct for winter solistice
  doy <- doy + 11

  # earths ecliptic
  j <- pi / 182.625
  axis <- (pi / 180) * 23.439

  # calculate daylength for all days
  dl <- lapply(doy, function(x){

    # Exposed radius part between sun's zenith and sun's circle
    m <- 1 - tan(latitude) * tan(axis * cos(j * x))

    # sun never appears or disappears
    if (m < 0) { m = 0 }
    if (m > 2) { m = 2 }

    # Exposed fraction of the sun's circle
    b <- acos(1 - m) / pi

    # Daylength (lat,day)
    return(b * 24)
  })

  # return the daylength vector
  return(unlist(dl))
}

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

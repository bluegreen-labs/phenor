#' Calculates AICc values
#'
#' Uses a set of measured and predicted values
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

rotate_cmip5 <- function(
  r = NULL,
  extent = c(-74, -65, 40, 48)
  ){

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
    m = raster::shift(raster::crop(r, extent2),-360)
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

#' Shapes model output
#'
#' The final output is based upon the class of the input data
#' and valid model estimates. Mainly, reshapes data to a spatial
#' raster format when required and traps NA and Inf values during
#' optimization schemes.
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
  if(methods::is(data, "phenor_map_data")){
    r = raster::raster(nrows = data$georeferencing$size[1],
                       ncols = data$georeferencing$size[2])
    raster::extent(r) <- data$georeferencing$extent
    sp::proj4string(r) <- sp::CRS(data$georeferencing$projection)
    r[] <- as.vector(doy)
    r[is.infinite(r)] <- NA
    return(r)
  } else {
    doy[is.na(doy)] <- 9999
    doy[is.infinite(doy)] <- 9999
    return(doy)
  }
}


#' Calculate daylength values
#'
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

#' Triangular temperature response function
#'
#' As defined in: Basler et al. 2016 (Agr. For. Meteorlogy)
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

#' Read parameter boundary values
#'
#' @param model a phenology model name
#' @param par_ranges location of the model parameter boundary file
#' @export
#' @examples
#'
#' T_response = triangular_temperature_response(T = 0:45)
#' \dontrun{
#' plot(0:45, T_response, type = "l")
#'}

pr_parameters <- function(
  model,
  par_ranges = system.file(
    "extdata",
    "parameter_ranges.csv",
    package = "phenor",
    mustWork = TRUE)
  ){

  if(missing(model)){
    stop("Please provide a model name")
  }

  # read in parameter ranges
  par_ranges = utils::read.table(par_ranges,
                                 header = TRUE,
                                 sep = ",")

  # subset the parameter range
  if (!any(par_ranges$model == model)){
    stop("parameters are not specified in the default parameter file.")
  }

  # extract parameter ranges is the model is available
  # in the file provided
  d = par_ranges[par_ranges$model == model,]

  d = d[,!is.na(d[1,])]
  d = d[,3:ncol(d)]
  d = as.matrix(d)

  # returns the converted temperature data
  return(d)
}

#' Outlier detection of transition dates
#'
#' Uses the 30-day rule (Schaber and Badeck, 2002) to remove outliers
#' from visual observations records. This function should be considered
#' for use with either PEP725 or USANPN data.
#'
#' @param df a nested list phenor formatted data file
#'
#' @return a nested list phenor formatted data file with outliers removed
#' @export

outlier_detection <- function(df){

  # check if this is a flat file or a nested one
  flat <- !is.null(df$site)

  if(flat){
    stop("
 This file is flattened after formating and can't be processed.
 If outlier detection is required, execute this step before
 combining nested list elements!
         ")
  }

  # Detect and remove outliers using the 30 day rule.
  # This routine mostly applies to visual observation data
  # in most other data streams outlier are removed with other routines.
  out <- lapply(df, function(x){

    # calculate the difference in days relative to
    # the true data
    diff_days <- abs(
      stats::predict(stats::lm(x$transition_dates ~ x$year)) -
        x$transition_dates
    )

    # flag all data which is further than 30 days
    # from the linear fit as suspect (outliers)
    # i.e. the 30 day rule
    x$transition_dates[which(diff_days >= 30)] <- NA

    s <- pr_fm_subset(x, !is.na(x$transition_dates))
    s$site <- x$site
    s$location <- x$location
    return(s)
  })

  # propagate class for post-processing
  class(out) = class(df)

  # return data
  return(out)
}

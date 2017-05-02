#' This file contains a model zoo of common phenological models.
#' All models take a data and parameter (par) input, but can vary
#' widely in terms of model structure. Models as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data: a nested list of data with on location:
#' 1. the date (doy or long format)
#' 2. the temperature data
#' 3. the photoperiod data (NA when not needed)
#' 4. a vector or matrix with necessary constants (NA when not needed)
#'    - long term mean temperature
#'    - latitude
#'    - etc...
#' @param par: a vector of parameter values, this is functions specific
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = TT(data = data, par = par)
#'}

# Basler, linear model (simple regression) FIX NOT CORRECT
LIN = function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  a = par[1]
  b = par[2]
  T_base = 5

  # convert DOY value to location in shifted
  # data set
  # t0 = which(data$doy == t0)

  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  tmean = data$tmean

  # convert vector to matrix of so provided
  if (is.vector(tmean)){
    tmean = as.matrix(tmean)
  }

  # calculate the model output on a column
  # by column basis (same for all other functions below)
  results = apply(tmean, 2, function(xt){

    # create forcing/chilling rate vector
    # forcing
    Rf = rep(0, length(xt))
    Rf[which(xt > T_base)] = xt[which(xt > T_base)] - T_base
    Rf[1:t0] = 0

    # sum all values from p0 (starting date)
    Sf = cumsum(Rf)

    # return the location(s) where the state of forcing
    # is exceeded by the critical threshold
    doy = data$doy[which(Sf >= F_crit)[1]]
    if (is.na(doy)){
      return(9999)
    } else {
      return(doy)
    }
  })
}

# Basler deepening rest model (NOT FINISHED CHECK)
DR = function(par, data, plot = TRUE){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 = par[1]
  a = par[2]
  b = par[3]
  c = par[4]
  T_base = par[5]

  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  tmean = data$tmean

  # convert vector to matrix of so provided
  if (is.vector(tmean)){
    tmean = as.matrix(tmean)
  }

  results = apply(tmean, 2, function(xt) {

    # grab length xt
    l = length(xt)

    # chilling
    Rc = rep(0, l)
    Rc[which(xt < T_opt & xt >= T_min)] = (xt[which(xt < T_opt & xt >= T_min)] - T_min)/
      (T_opt - T_min)
    Rc[which(xt < T_max & xt >= T_opt)] = (xt[which(xt < T_max & xt >= T_opt)] - T_opt)/
      (T_max - T_opt)
    Rc[1:t0] = 0
    Sc = cumsum(Rc)

    # k
    k = rep(0, l)
    k[which(Sc >= C_req)] = 1

    # forcing
    Rf = rep(0, l)
    Rf[which(xt > T_base)] = xt[which(xt > T_base)] - T_base
    Rf = Rf * k
    Rf[1:t0] = 0
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(data$doy[which(Sf >= (w * exp( f * Sc)) )[1]])
  })
}

DRb = function(par, data, plot = TRUE){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 = par[1]
  a = par[2]
  b = par[3]
  c = par[4]
  T_base = par[5]

  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  tmean = data$tmean

  # convert vector to matrix of so provided
  if (is.vector(tmean)){
    tmean = as.matrix(tmean)
  }

  results = apply(tmean, 2, function(xt) {

    # grab length xt
    l = length(xt)

    # chilling
    Rc = rep(0, l)
    Rc[which(xt < T_opt & xt >= T_min)] = (xt[which(xt < T_opt & xt >= T_min)] - T_min)/
      (T_opt - T_min)
    Rc[which(xt < T_max & xt >= T_opt)] = (xt[which(xt < T_max & xt >= T_opt)] - T_opt)/
      (T_max - T_opt)
    Rc[1:t0] = 0
    Sc = cumsum(Rc)

    # k
    k = rep(0, l)
    k[which(Sc >= C_req)] = 1

    # forcing
    Rf = rep(0, l)
    Rf[which(xt > T_base)] = xt[which(xt > T_base)] - T_base
    Rf = Rf * k
    Rf[1:t0] = 0
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(data$doy[which(Sf >= (w * exp( f * Sc)) )[1]])
  })
}

# Basler Four Phase Model
#FP
#FPb

# Basler DORMPHOT model
DR = function(par, data, plot = TRUE){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 = par[1]
  a = par[2]
  b = par[3]
  c = par[4]
  T_base = par[5]

  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  tmean = data$tmean

  # convert vector to matrix of so provided
  if (is.vector(tmean)){
    tmean = as.matrix(tmean)
  }

  results = apply(tmean, 2, function(xt) {

    # grab length xt
    l = length(xt)

    # chilling
    Rc = rep(0, l)
    Rc[which(xt < T_opt & xt >= T_min)] = (xt[which(xt < T_opt & xt >= T_min)] - T_min)/
      (T_opt - T_min)
    Rc[which(xt < T_max & xt >= T_opt)] = (xt[which(xt < T_max & xt >= T_opt)] - T_opt)/
      (T_max - T_opt)
    Rc[1:t0] = 0
    Sc = cumsum(Rc)

    # k
    k = rep(0, l)
    k[which(Sc >= C_req)] = 1

    # forcing
    Rf = rep(0, l)
    Rf[which(xt > T_base)] = xt[which(xt > T_base)] - T_base
    Rf = Rf * k
    Rf[1:t0] = 0
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(data$doy[which(Sf >= (w * exp( f * Sc)) )[1]])
  })
}


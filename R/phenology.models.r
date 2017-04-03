#' This file contains a model zoo of common phenological models.
#' All models take a data and parameter (par) input, but can vary
#' widely in terms of model structure. Models as defined in
#' Melaas et al. 2016 (Global Change Biology)
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
#' estimate <- SEQ1.3(data,par)
#'
#' # estimate will be an estimated timing of a phenophase
#' # the format depends on the date as forwarded into the
#' # model
#'}

#--- Spring warming models

SW1 = function(par, data, plot = FALSE){

  # exit the routine as some parameters are missing
  if (length(par) != 2){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 = par[1]
  F_star = par[2]

  # split out the temperature data
  date = data$doy

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

    # grab length xt
    l = length(xt)

    # exponential conversion
    xt_exp = 28.4 / (1 + exp(3.4 - 0.185 * xt))

    # create forcing/chilling rate vector
    # and fill it with the necessary values
    # on locations fullfiling the forcing/chilling requirement
    Rf = rep(0, l)
    Rf[which(xt > 0)] = xt_exp[which(xt > 0)]

    # remove all values before t0
    Rf[1:t0] = 0

    # sum all values from p0 (starting date)
    Sf = cumsum(Rf)

    # return the location(s) where the state of forcing
    # is exceeded by the critical threshold
    return(date[which(Sf >= F_star)[1]])
  })
}

SW2 = function(par, data, plot = FALSE){

  # exit the routine as some parameters are missing
  if (length(par) != 2){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  p0 = par[1]
  F_star = par[2]

  # split out the temperature data
  date = data$doy

  # latitude
  lat = data$location[1]

  # doy
  doy = data$doy

  # data
  tmean = data$tmean

  # calculation of the daylight hours
  photoperiod = unlist(daylength(doy,lat)[1])

  # convert vector to matrix of so provided
  if (is.vector(tmean)){
    tmean = as.matrix(tmean)
  }

  results = apply(tmean, 2, function(xt) {

    # grab length xt
    l = length(xt)

    # exponential conversion
    xt_exp = 28.4 / (1 + exp(3.4 - 0.185 * xt))

    # create forcing/chilling rate vector
    # and fill it with the necessary values
    # on locations fullfiling the forcing/chilling requirement
    Rf = rep(0, l)
    Rf[which(xt > 0)] = xt_exp[which(xt > 0)]

    # find where the photoperiod threshold is
    # exceeded
    p0_loc = which(photoperiod > p0)

    if (length(p0_loc) == 0){
      return(9999)
    }

    Rf[1:p0_loc] = 0

    # sum all values from p0 (starting date)
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(date[which(Sf >= F_star)[1]])
  })
}

SW3 = function(par, data, plot = FALSE){

  # exit the routine as some parameters are missing
  if (length(par) != 3){
    print(par)
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 = par[1]
  a = par[2]
  b = par[3]

  # split out the temperature data
  doy = data$doy

  # grab data matrix
  tmean = data$tmean

  # check formulation
  Tbar = mean(data$ltm, na.rm = TRUE)

  # convert vector to matrix of so provided
  if (is.vector(tmean)){
    tmean = as.matrix(tmean)
  }

  results = apply(tmean, 2, function(xt) {
    # grab length xt
    l = length(xt)

    # exponential conversion
    xt_exp = 28.4 / (1 + exp(3.4 - 0.185 * xt))

    # create forcing/chilling rate vector
    # and fill it with the necessary values
    # on locations fullfiling the forcing/chilling requirement
    Rf = rep(0, l)
    Rf[which(xt > 0)] = xt_exp[which(xt > 0)]
    Rf[1:t0] = 0

    # sum all values from p0 (starting date)
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(doy[which(Sf >= (a * Tbar + b))[1]])
  })
}

SW4 = function(par, data, plot = FALSE){

  # exit the routine as some parameters are missing
  if (length(par) < 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  p0 = par[1]
  a = par[2]
  b = par[3]
  c = par[4]
  d = par[5]

  # calculation of the daylight hours
  photoperiod = unlist(daylength(data$doy,
                                 data$location[1])[1])

  # mean annual temperature
  Tbar = mean(data$ltm)

  results = apply(data$tmean, 2, function(xt) {
    # grab length xt
    l = length(xt)

    # exponential conversion
    xt_exp = 28.4 / (1 + exp(3.4 - 0.185 * xt))

    # create forcing/chilling rate vector
    # and fill it with the necessary values
    # on locations fullfiling the forcing/chilling requirement
    Rf = rep(0, l)
    Rf[which(xt > 0)] = xt_exp[which(xt > 0)]

    # find where the photoperiod threshold is
    # exceeded
    p0_loc = which(photoperiod > p0)
    if(length(p0_loc) == 0){
      return(9999)
    } else {
      Rf[1:p0_loc] = 0
    }

    # sum all values from p0 (starting date)
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(data$doy[which(Sf >= (a * Tbar + b))[1]])
  })
}

#---Sequential models

SEQ1 = function(par, data, plot = TRUE){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 = par[1]
  F_star = par[2]
  C_star = par[3]
  Tf = par[4]
  Tc = par[5]

  results = apply(data$tmean, 2, function(xt) {

    # grab length xt
    l = length(xt)

    # chilling
    Rc = rep(0, l)
    Rc[which(xt < Tc)] = 1
    Rc[1:t0] = 0
    Sc = cumsum(Rc)
    t1 = which(Sc >= C_star)[1]

    # trap NA value for t1
    if(is.na(t1) | length(t1) == 0 ){
      return(9999)
    }

    # forcing
    Rf = rep(0, l)
    Rf[which(xt > Tf)] = xt[which(xt > Tf)] - Tf
    Rf[1:t1] = 0
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(data$doy[which(Sf >= F_star)[1]])
  })
}

SEQ2 = function(par, data, plot = TRUE){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  p0 = par[1]
  F_star = par[2]
  C_star = par[3]
  Tf = par[4]
  Tc = par[5]

  # calculation of the daylight hours
  photoperiod = unlist(daylength(data$doy,data$location[1])[1])

  results = apply(data$tmean, 2, function(xt) {

    # grab length xt
    l = length(xt)

    # chilling
    Rc = rep(0, l)
    Rc[which(xt < Tc)] = 1

    # find where the photoperiod threshold is
    # exceeded
    p0_loc = which(photoperiod < p0)
    if(length(p0_loc) == 0 ){
      return(9999)
    }
    Rc[!p0_loc] = 0

    Sc = cumsum(Rc)
    t1 = which(Sc >= C_star)[1]

    if(is.na(t1) | length(t1) == 0 | p0_loc[1] > t1 ){
      return(9999)
    }

    # forcing
    Rf = rep(0, l)
    Rf[which(xt > Tf)] = xt[which(xt > Tf)] - Tf
    Rf[1:t1] = 0
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(data$doy[which(Sf >= F_star)[1]])
  })
}

SEQ3 = function(data = data, par = par, plot = FALSE){

  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  p0 = par[1]
  C_star = par[2]
  Tf = par[3]
  Tc = par[4]
  a = par[5]
  b = par[6]

  # calculation of the daylight hours
  photoperiod = unlist(daylength(data$doy,data$location[1])[1])

  # mean annual temperature
  Tbar = mean(data$ltm)

  results = apply(data$tmean, 2, function(xt) {

    # grab length xt
    l = length(xt)

    # chilling
    Rc = rep(0, l)
    Rc[which(xt < Tc)] = 1

    # find where the photoperiod threshold is
    # exceeded
    p0_loc = which(photoperiod < p0)
    if(length(p0_loc) == 0 ){
      return(9999)
    }
    Rc[!p0_loc] = 0

    Sc = cumsum(Rc)
    t1 = which(Sc >= C_star)[1]

    if(is.na(t1)){
      return(9999)
    }

    # forcing
    Rf = rep(0, l)
    Rf[which(xt > Tf)] = xt[which(xt > Tf)] - Tf
    Rf[1:t1] = 0
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(data$doy[which(Sf >= (a*Tbar + b))[1]])
  })
}

SEQ4 = function(par, data, plot = FALSE) {

  # exit the routine as some parameters are missing
  if (length(par) != 8){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  p0 = par[1]
  C_star = par[2]
  Tf = par[3]
  Tc = par[4]
  a = par[5]
  b = par[6]
  c = par[7]
  d = par[8]

  # calculation of the daylight hours
  photoperiod = unlist(daylength(data$doy,data$location[1])[1])

  # mean annual temperature
  Tbar = mean(data$ltm)

  results = apply(data$tmean, 2, function(xt) {

    # grab length xt
    l = length(xt)

    # chilling
    Rc = rep(0, l)
    Rc[which(xt < Tc)] = 1

    # find where the photoperiod threshold is
    # exceeded
    p0_loc = which(photoperiod < p0)
    if(length(p0_loc) == 0 ){
      return(9999)
    }
    Rc[!p0_loc] = 0

    Sc = cumsum(Rc)
    t1 = which(Sc >= a*Tbar + b)[1]

    if(is.na(t1)){
      return(9999)
    }

    # forcing
    Rf = rep(0, l)
    Rf[which(xt > Tf)] = xt[which(xt > Tf)] - Tf
    Rf[1:t1] = 0
    Sf = cumsum(Rf)

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(data$doy[which(Sf >= (c*Tbar + d))[1]])
  })
}

#--- frost hardiness models

# frost hardiness model, needs F* parameter to set threshold
# for phenology first / second order model of frost-hardiness
# fluctuating development form Leinonen et al. 1995)
# Model D in A framework for modelling the annual
# cycle of trees in boreal and temperate regions

C1 = function(par,data){

  # average temperature
  T = data$tmean

  # length T
  l = length(T)

  # get latitude
  lat = data$location[1]

  # split out the temperature data
  doy = data$doy

  # calculate length of the night
  NL = 24 - unlist(daylength(doy, lat)[1])

  # split out parameters into
  # readable format
  a = par[1]
  b = par[2]
  t = par[3]
  t2 = par[4]
  T8 = 11.3 # parameter
  Rmin = -5 # parameter

  # create empty state vectors
  Ssh = rep(0,l)
  Rah = rep(0,l)
  Rh = rep(0,l)
  Sh = rep(0,l)

  NL1 = 10
  NL2 = 15
  deltaRtMax = -47
  deltaRpMax = -18.5
  Rmin = -4.5

  # calculate the stationary frost hardiness (D1)
  # based upon temperature
  Rt = ifelse(T <= T8, a * T + b, a * T8 + b)
  Rp = ifelse(NL1 <= NL & NL <= NL2,
              NA,
              ifelse(NL < NL1,
                     0,
                     deltaRpMax))

  if (all(is.na(Rp))){
    Rp = 0
  }else{
    Rp = zoo::na.approx(Rp)
  }

  # sum driving factors
  Ssh = Rmin + Rt + Rp

  # calculate the rate of change of frost hardiness
  # and asymptotic frost hardiness (C2/3 and D2/3,
  # where D3 equals the sum of Rah/Rh up until t/i).
  for (i in 1:l){
    Sh = sum(Rh[1:i])
    Sah = sum(Rah[1:i])
    Rh[i] = 1/t * (Ssh[i] - Sh)
    Rah[i] = 1/t2 * (Ssh[i] - Sah)
  }
  return(cumsum(Rh))
}

#Second order model of frost hardiness
D1 = function(par,data){

  # average temperature
  T = data$tmean

  # length T
  l = length(T)

  # get latitude
  lat = data$location[1]

  # split out the temperature data
  doy = data$doy

  # calculate length of the night
  NL = 24 - unlist(daylength(doy, lat)[1])

  # split out parameters into
  # readable format
  a = par[1]
  b = par[2]
  t = par[3]
  t2 = par[4]
  T8 = 11.3 # parameter
  Rmin = -5 # parameter

  # create empty state vectors
  Ssh = rep(0,l)
  Rah = rep(0,l)
  Rh = rep(0,l)
  Sh = rep(0,l)

  NL1 = 10
  NL2 = 15
  deltaRtMax = -47
  deltaRpMax = -18.5
  Rmin = -4.5

  # calculate the stationary frost hardiness (D1)
  # based upon temperature
  Rt = ifelse(T <= T8, a * T + b, a * T8 + b)
  Rp = ifelse(NL1 <= NL & NL <= NL2,
              NA,
              ifelse(NL < NL1,
                     0,
                     deltaRpMax))

  if (all(is.na(Rp))){
    Rp = 0
  }else{
    Rp = zoo::na.approx(Rp)
  }

  # sum driving factors
  Ssh = Rmin + Rt + Rp

  # calculate the rate of change of frost hardiness
  # and asymptotic frost hardiness (C2/3 and D2/3,
  # where D3 equals the sum of Rah/Rh up until t/i).
  for (i in 1:l){
    Sh = sum(Rh[1:i])
    Sah = sum(Rah[1:i])
    Rh[i] = 1/t * (Ssh[i] - Sh)
    Rah[i] = 1/t2 * (Ssh[i] - Sah)
  }
  return(cumsum(Rh))
}

#--- grassland pulse response model

# simple pollen model, good substitute in order to include precip
# Garcia-Mozo et al. 2009 AFM
# requires time series which do not start in the previous year
# and run form Jan 1 to Dec 31
GR1 = function(data,
               par = c(-1.92, -3.94, 70.59, 9.6, 9.9),
               plot = FALSE){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  p0 = par[1]
  F_star = par[2]
  R_star = par[3]
  a = par[4]
  b = par[5]

  # split out the temperature data
  tmean = data$tmean

  # get precip data
  precip = data$precip

  # split out the temperature data
  doy = data$doy

  # get latitude
  lat = data$location[1]

  # calculation of the daylight hours
  photoperiod = unlist(daylength(doy,lat)[1])

  # create temporary matrix
  tmp_data = rbind(tmean,precip)

  results = apply(tmp_data,2,function(x){

    # subset the data
    xt = x[1:365]

    prcp = x[366:length(x)]

    # format variables
    l = length(xt)
    Rt = rep(0, l)

    # which day sets the accumulation of
    # precipitation, technically the doy
    # subset should not be used, but for
    # clarity I do so anyway (same below)
    t0 = doy[which(photoperiod > p0)[1]]

    # trap errors
    if (is.null(t0) | is.na(t0)){
      return(9999)
    }

    # for locations Rp, set precip to 0 which
    # blocks further development of the model
    prcp[doy <= t0] = 0

    # now apply a rolling sum over the precip data
    Rt = zoo::rollapply(
      data = prcp,
      width = 7,
      FUN = sum,
      align = "left",
      fill = NA
    )

    # find doy from which to accumulate temperature
    t1 = doy[which(Rt > R_star)[1]]

    # set all temperatures before t1 to 0
    # no accumulation
    xt[doy < t1] = 0

    # takes daily temperature values as input
    # calculate the difference with the base temperature
    Rf = 1 / (1 + exp(a * (xt + b)))

    # sum all values from p0 (starting date)
    Sf = cumsum(Rf)

    if (plot){
      par(mfrow=c(2,1))
      plot(Sf)
      lines(xt)
      abline(h = F_star)
      abline(v = t0)
      plot(prcp)
    }

    # return the location where the state of forcing
    # is exceeded by the critical threshold
    return(doy[which(Sf >= F_star )[1]])
  })
  return(results)
}

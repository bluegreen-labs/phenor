#' Calculates estimated phenophases for a set of common phenology
#' models as listed in Melaas et al. 2013. In particular it uses a
#' photoperiod constrain to account for changes in latitude instead
#' of explicitely changing parameter p0 (see paper).
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
#' estimate <- SEQ1.3(data,par)
#' 
#' # estimate will be an estimated timing of a phenophase
#' # the format depends on the date as forwarded into the
#' # model
#' 

# Spring warming models
SW1.1 = function(par,data){ 
  
  # exit the routine as some parameters are missing
  if (length(par) < 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument
  Tf = par[1]
  p0 = par[2]
  F = par[3]

  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  
  # takes daily temperature values as input
  # calculate the difference with the base temperature
  xt = xt - Tf
  
  # grab length xt
  l = length(xt)
  
  # create forcing/chilling rate vector
  # and fill it with the necessary values
  # on locations fullfiling the forcing/chilling requirement
  Rf = rep(0,l)
  Rf[which(xt > Tf)] = xt[which(xt > Tf)] 
  
  # sum all values from p0 (starting date)
  Sf = cumsum(Rf[p0:l])
  
  # return the location where the state of forcing
  # is exceeded by the critical threshold
  return(date[which[Sf >= F][1]])
}

SW1.2 = function(par,data){ 
  
  # exit the routine as some parameters are missing
  if (length(par) < 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument
  Tf = par[1]
  p0 = par[2]
  F = par[3]
  a = par[4]
  b = par[5]
  Tbar = par[6]
  
  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  
  # takes daily temperature values as input
  # calculate the difference with the base temperature
  xt = xt - (a*Tbar +b)
  
  # grab length xt
  l = length(xt)
  
  # create forcing/chilling rate vector
  # and fill it with the necessary values
  # on locations fullfiling the forcing/chilling requirement
  Rf = rep(0,l)
  Rf[which(xt > (a*Tbar + b))] = xt[which(xt > (a*Tbar + b))] 
  
  # sum all values from p0 (starting date)
  Sf = cumsum(Rf[p0:l])
  
  # return the location where the state of forcing
  # is exceeded by the critical threshold
  return(date[which[Sf >= F][1]])
}

SW1.3 = function(par,data){ 
  
  # exit the routine as some parameters are missing
  if (length(par) < 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument
  Tf = par[1]
  p0 = par[2]
  F = par[3]
  a = par[4]
  b = par[5]
  Tbar = par[6]
  
  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  
  # takes daily temperature values as input
  # calculate the difference with the base temperature
  xt = xt - Tf
  
  # grab length xt
  l = length(xt)
  
  # create forcing/chilling rate vector
  # and fill it with the necessary values
  # on locations fullfiling the forcing/chilling requirement
  Rf = rep(0,l)
  Rf[which(xt > Tf)] = xt[which(xt > Tf)] 
  
  # sum all values from p0 (starting date)
  Sf = cumsum(Rf[p0:l])
  
  # return the location where the state of forcing
  # is exceeded by the critical threshold
  return(date[which[Sf >= (a*Tbar + b)][1]])
}

SW1.4 = function(par,data){ 
  
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
  Tbar = par[6]
  
  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  
  # takes daily temperature values as input
  # calculate the difference with the base temperature
  xt = xt - (a*Tbar + b)
  
  # grab length xt
  l = length(xt)
  
  # create forcing/chilling rate vector
  # and fill it with the necessary values
  # on locations fullfiling the forcing/chilling requirement
  Rf = rep(0,l)
  Rf[which(xt > (a*Tbar + b))] = xt[which(xt > (a*Tbar + b))] 
  
  # sum all values from p0 (starting date)
  Sf = cumsum(Rf[p0:l])
  
  # return the location where the state of forcing
  # is exceeded by the critical threshold
  return(date[which[Sf >= (c*Tbar + d)][1]])
}

SW2.1 = function(par,data){ 
  
  # exit the routine as some parameters are missing
  if (length(par) < 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument
  p0 = par[1]
  F = par[2]
  
  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  
  # takes daily temperature values as input
  # calculate the difference with the base temperature
  Rf = 28.4/(1-exp(-0.185*xt-18.4))
  
  # grab length xt
  l = length(xt)
  
  # create forcing/chilling rate vector
  # and fill it with the necessary values
  # on locations fullfiling the forcing/chilling requirement
  Rf[which(xt > 0)] = xt[which(xt > 0)]
  Rf[which(xt <= 0)] = 0
  
  # sum all values from p0 (starting date)
  Sf = cumsum(Rf[p0:l])
  
  # return the location where the state of forcing
  # is exceeded by the critical threshold
  return(date[which[Sf >= (c*Tbar + d)][1]])
}

SW2.2 = function(par,data){ 
  
  # exit the routine as some parameters are missing
  if (length(par) < 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument
  p0 = par[1]
  a = par[3]
  b = par[4]
  Tbar = par[5]
  
  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  
  # takes daily temperature values as input
  # calculate the difference with the base temperature
  Rf = 28.4/(1-exp(-0.185*xt-18.4))
  
  # create forcing/chilling rate vector
  # and fill it with the necessary values
  # on locations fullfiling the forcing/chilling requirement
  Rf[which(xt > 0)] = xt[which(xt > 0)] 
  Rf[which(xt <= 0)] = 0
  
  # sum all values from p0 (starting date)
  Sf = cumsum(Rf[p0:l])
  
  # return the location where the state of forcing
  # is exceeded by the critical threshold
  return(date[which[Sf >= (a*Tbar + b)][1]])
}

# Sequential models
SEQ1.1 = function(data=data,par=par){
  
  # exit the routine as some parameters are missing
  if (length(par) < 7){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument
  Tc = par[1]
  Tf = par[2]
  C = par[3]
  F = par[4]
  p0 = par[5]
  t1 = par[6]
  
  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  photoperiod = data[[3]]
  
  # put constants in the data frame,
  # this includes latitude, mean annual temperature etc
  # one row per site, similar to the temperature data
  Tbar = data[[4]]
  
  # grab length xt
  l = length(xt)
  
  Rc = rep(0,l)
  Rc[which(xt < Tc)] = 1
  
  # subtract the threshold temperature
  xt = xt - Tf # calculate the difference with the threshold
  Rf = xt # assign a new variable name
  Rf[Rf<0] = 0 # set all negative values to 0
  
  # Integrations over time for both Rc and Rf:
  # accumulation of chilling degree days with the range defined
  # by the date where p0 == the location's photoperiod / daylength
  # (in daylength hours)

  # photoperiod requirement of integration time
  # of the chilling degree days
  if (min(photoperiod) >= p0){
    loc = 113
  }else{
    if(max(photoperiod) <= p0){
      loc = 1
    }else{
      loc = which(photoperiod <= p0)[1]
    }
  }
  
  # set values before p0 to 0 (don't count in cumsum())
  Rc[1:loc] = 0
  Sc = cumsum(Rc)
  
  # starting date for accumulation of forcing
  t1 = which(Sc >= C)[1] # pick the first one that exceeds the C threshold hence [1]
  
  # trap NA value for t1
  if (is.na(t1)){
    return(NA)
  }
  
  # set values before t1 to 0 (don't count in cumsum())
  Rf[1:t1] = 0
  Sf = cumsum(Rf)
  
  # find your phenological metric
  phenological_metric = date[which(Sf >= F)[1]]
  
  if(is.na(phenological_metric)){
    return(NA)
  }else{
    return(phenological_metric)
  }
}


SEQ1.2 = function(data=data,par=par){
  
  # exit the routine as some parameters are missing
  if (length(par) < 7){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument
  C = par[3]
  F = par[4]
  p0 = par[5]
  a = par
  b = par
  c = par
  d = par
  Tbar = par
  
  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  photoperiod = data[[3]]
  
  # put constants in the data frame,
  # this includes latitude, mean annual temperature etc
  # one row per site, similar to the temperature data
  Tbar = data[[4]]
  
  # grab length xt
  l = length(xt)
  
  Rc = rep(0,l)
  Rc[which(xt < (c*Tbar + d))] = 1
  
  # subtract the threshold temperature
  xt = xt - (a*Tbar + b) # calculate the difference with the threshold
  Rf = xt # assign a new variable name
  Rf[Rf<0] = 0 # set all negative values to 0
  
  # Integrations over time for both Rc and Rf:
  # accumulation of chilling degree days with the range defined
  # by the date where p0 == the location's photoperiod / daylength
  # (in daylength hours)

  # photoperiod requirement of integration time
  # of the chilling degree days
  if (min(photoperiod) >= p0){
    loc = 113
  }else{
    if(max(photoperiod) <= p0){
      loc = 1
    }else{
      loc = which(photoperiod <= p0)[1]
    }
  }
  
  # set values before p0 to 0 (don't count in cumsum())
  Rc[1:loc] = 0
  Sc = cumsum(Rc)
  
  # starting date for accumulation of forcing
  t1 = which(Sc >= C)[1] # pick the first one that exceeds the C threshold hence [1]
  
  # trap NA value for t1
  if (is.na(t1)){
    return(NA)
  }
  
  # set values before t1 to 0 (don't count in cumsum())
  Rf[1:t1] = 0
  Sf = cumsum(Rf)
  
  # find your phenological metric
  phenological_metric = date[which(Sf >= F)[1]]
  
  if(is.na(phenological_metric)){
    return(NA)
  }else{
    return(phenological_metric)
  }
}

SEQ1.3 = function(data=data,par=par){
  
  # exit the routine as some parameters are missing
  if (length(par) < 7){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument
  Tc = par[1]
  Tf = par[2]
  C = par[3]
  p0 = par[4]
  a = par[5]
  b = par[6]
  tp = par[8]
  
  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  photoperiod = data[[3]]
  
  # put constants in the data frame,
  # this includes latitude, mean annual temperature etc
  # one row per site, similar to the temperature data
  Tbar = data[[4]]
  
  # grab length xt
  l = length(xt)
  
  Rc = rep(0,l)
  Rc[which(xt < Tc)] = 1
  
  # subtract the threshold temperature
  xt = xt - Tf # calculate the difference with the threshold
  Rf = xt # assign a new variable name
  Rf[Rf<0] = 0 # set all negative values to 0
  
  # Integrations over time for both Rc and Rf:
  # accumulation of chilling degree days with the range defined
  # by the date where p0 == the location's photoperiod / daylength
  # (in daylength hours)

  # photoperiod requirement of integration time
  # of the chilling degree days
  if (min(photoperiod) >= p0){
    loc = 113
  }else{
    if(max(photoperiod) <= p0){
      loc = 1
    }else{
      loc = which(photoperiod <= p0)[1]
    }
  }
  
  # set values before p0 to 0 (don't count in cumsum())
  Rc[1:loc] = 0
  Sc = cumsum(Rc)
  
  # starting date for accumulation of forcing
  t1 = which(Sc >= C)[1] # pick the first one that exceeds the C threshold hence [1]
  
  if (is.na(t1)){
    #t1 = round(tp)
    return(NA)
  }
  
  # set values before t1 to 0 (don't count in cumsum())
  Rf[1:t1] = 0
  Sf = cumsum(Rf)
  
  # find your phenological metric
  phenological_metric = date[which(Sf >= (a*Tbar + b))[1]]
  
  if(is.na(phenological_metric)){
    return(0)
  }else{
    return(phenological_metric)
  }
}

SEQ1.4 = function(data=data,par=par){
  
  # exit the routine as some parameters are missing
  if (length(par) < 7){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument
  Tc = par[1]
  Tf = par[2]
  F = par
  C = par[3]
  p0 = par[4]
  a = par[5]
  b = par[6]
  tp = par[8]
  
  # split out the temperature data
  date = data[[1]]
  
  # matrix nr_obs x yearly values between 21 sept. previous year
  # and 21 sept. current year
  xt = data[[2]]
  photoperiod = data[[3]]
  
  # put constants in the data frame,
  # this includes latitude, mean annual temperature etc
  # one row per site, similar to the temperature data
  Tbar = data[[4]]
  
  # grab length xt
  l = length(xt)
  
  Rc = rep(0,l)
  Rc[which(xt < Tc)] = 1
  
  # subtract the threshold temperature
  xt = xt - Tf # calculate the difference with the threshold
  Rf = xt # assign a new variable name
  Rf[Rf<0] = 0 # set all negative values to 0
  
  # Integrations over time for both Rc and Rf:
  # accumulation of chilling degree days with the range defined
  # by the date where p0 == the location's photoperiod / daylength
  # (in daylength hours)
  
  # photoperiod requirement of integration time
  # of the chilling degree days
  if (min(photoperiod) >= p0){
    loc = 113
  }else{
    if(max(photoperiod) <= p0){
      loc = 1
    }else{
      loc = which(photoperiod <= p0)[1]
    }
  }
  
  # set values before p0 to 0 (don't count in cumsum())
  Rc[1:loc] = 0
  Sc = cumsum(Rc)
  
  # starting date for accumulation of forcing
  t1 = which(Sc >= C)[1] # pick the first one that exceeds the C threshold hence [1]
  
  if (is.na(t1)){
    #t1 = round(tp)
    return(NA)
  }
  
  # set values before t1 to 0 (don't count in cumsum())
  Rf[1:t1] = 0
  Sf = cumsum(Rf)
  
  # find your phenological metric
  phenological_metric = date[which(Sf >= (a*Tbar + b))[1]]
  
  if(is.na(phenological_metric)){
    return(0)
  }else{
    return(phenological_metric)
  }
}
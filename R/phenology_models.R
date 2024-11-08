#' Alternating model as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = AT(data = data, par = par)
#'}

AT <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  a <- par[3]
  b <- par[4]
  c <- par[5]

  # chilling
  Rc <- ifelse(data$Ti >= T_base, 0, 1)
  Rc[1:t0,] <- 0
  Sc <- apply(Rc, 2, cumsum)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf <= 0] <- 0
  Rf[1:t0,] <- 0
  Sf <- apply(Rf, 2, cumsum)

  Sfc <- Sf - (a + b * exp(c * Sc))

  doy <- apply(Sfc, 2, function(x){
    data$doy[which(x > 0)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Chilling Degree Day (CDD) adapted from
#' Jeong & Medvigny 2014 (Global Ecology & Biogeography)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- CDD(data = data, par = par)
#'}

CDD <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]

  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' DormPhot model as defined in
#' Caffarra, Donnelly and Chuine 2011 (Clim. Res.)
#' parameter ranges are taken from Basler et al. 2016
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = DP(data = data, par = par)
#'}

DP <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 11){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  a <- par[1]
  b <- par[2]
  c <- par[3]
  d <- par[4]
  e <- par[5]
  f <- par[6]
  g <- par[7]
  F_crit <- par[8]
  C_crit <- par[9]
  L_crit <- par[10]
  D_crit <- par[11]

  # set the t0 value if necessary (~Sept. 1)
  t0 <- which(data$doy < 1 & data$doy == -121)

  # dormancy induction
  # (this is all vectorized doing cell by cell multiplications with
  # the sub matrices in the nested list)
  DR <- 1/(1 + exp(a * (data$Ti - b))) * 1/(1 + exp(10 * (data$Li - L_crit)))
  if (!length(t0) == 0){
    DR[1:t0,] <- 0
  }
  DS <- apply(DR,2, cumsum)

  # chilling
  CR <- 1/(1 + exp(c * (data$Ti - d)^2 + (data$Ti - d) ))
  CR[DS < D_crit] <- 0
  CS <- apply(CR,2, cumsum)

  # forcing
  dl50 <- 24 / (1 + exp(f * (CS - C_crit))) # f >= 0
  t50 <- 60 / (1 + exp(g * (data$Li - dl50))) # g >= 0
  Rf <- 1/(1 + exp(e * (data$Ti - t50))) # e <= 0
  Rf[DS < D_crit] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' DU drought model
#'
#' Floral induction model by
#' Chen et al. 2017 (Journal of Ecology)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model, tropical forest
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- DU(data = data, par = par)
#'}

DU <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  ni <- round(par[1])
  nd <- round(par[2])
  P_base <- par[3]
  F_crit <- par[4]

  # calculate floral induction time series
  Fd <- apply(data$Pi, 2, function(x){
    Fi <- zoo::rollapply(x, ni, mean, align = "right", fill = 0, na.rm = TRUE)
    Fd <- c(rep(0,nd),Fi[1:(length(Fi) - nd)])
  })

  # threshold value
  Fd <- Fd - P_base
  Fd[Fd < 0] <- 0

  # DOY of budburst criterium
  doy <- apply(Fd, 2, function(xt){
    data$doy[which(xt >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Thermal Time grassland pollen model
#'
#' Pulse response model on precipitation as defined in
#' Garcia-Mozo et al. 2009 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- GRP(data = data, par = par)
#'}

GRP <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  b <- par[1]
  c <- par[2]
  F_crit <- par[3]
  P_crit <- par[4]
  L_crit <- par[5]

  # Light requirement must be met
  # and daylength increasing (as per original specs)
  P_star <- ifelse(data$Li >= L_crit, 1, 0)
  P_star[diff(data$Li) < 0] <- 0

  # rainfall accumulation
  data$Pi <- data$Pi * P_star
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)

  # This avoids inefficient moving window
  # approaches which are computationally
  # expensive (break the routine when the
  # criterium is met instead of a full run
  # for a year)

  # assign output matrix
  R_star <- matrix(0,rows,cols)

  # fill in values where required
  # until P_crit is met
  for (i in 1:cols){
    for (j in 1:rows){
      if (j == rows - 7){
        break
      }
      if(sum(data$Pi[j:(j+7),i]) >= P_crit){
        R_star[j:rows,i] <- 1
        break
      }
    }
  }

  # create forcing/chilling rate vector
  # forcing
  Rf <- 1 / (1 + exp(b * (data$Ti + c))) * R_star

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Linear model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @param spring a vector defining spring as a couple of DOY values
#' default is March, April, May or DOY 60 - 151
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- LIN(data = data, par = par)
#'}

LIN <- function(par, data, spring = c(60,151)){

  # exit the routine as some parameters are missing
  if (length(par) != 2){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  a <- par[1]
  b <- par[2]

  # calculate location of "spring" values
  spring_loc <- data$doy %in% seq(spring[1],spring[2])

  # calculate the mean temperature for this range
  mean_spring_temp <- apply(data$Ti,2,function(xt){
    mean(xt[spring_loc])
  })

  # linear regression
  doy <- a * mean_spring_temp + b

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' M1 model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- M1(data = data, par = par)
#'}

M1 <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- par[1]
  T_base <- par[2]
  k <- par[3]
  F_crit <- par[4]

  # create forcing/chilling rate vector
  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- ((data$Li / 10) ^ k) * Rf
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' M1s model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' with a sigmoidal temperature response (Kramer 1994)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- M1s(data = data, par = par)
#'}

M1s <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- par[1]
  b <- par[2]
  c <- par[3]
  k <- par[4]
  F_crit <- par[5]

  # create forcing/chilling rate vector
  # forcing
  Rf <- 1 / (1 + exp(-b * (data$Ti - c)))
  Rf <- ((data$Li / 10) ^ k) * Rf
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Null model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' returns the mean across all validation dates
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- null(data = data)
#'}

null <- function(data){
  rep(round(mean(data$transition_dates,na.rm=TRUE)),
      length(data$transition_dates))
}

#' Parallel model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = PA(data = data, par = par)
#'}

PA <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 9){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 <- round(par[1])
  t0_chill <- round(par[2])
  T_base <- par[3]
  T_opt <- par[4]
  T_min <- par[5]
  T_max <- par[6]
  C_ini <- par[7]
  F_crit <- par[8]
  C_req <- par[9]

  # chilling
  Rc <- triangular_temperature_response(data$Ti,
                                       T_opt = T_opt,
                                       T_min = T_min,
                                       T_max = T_max)
  Rc[1:t0_chill,] <- 0
  Sc <- apply(Rc, 2, cumsum)

  k <- C_ini + Sc * (1 - C_ini)/C_req
  k[Sc >= C_req] = 1

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Parallel model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' using a bell shaped chilling response
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- PAb(data = data, par = par)
#'}

PAb <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 9){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  t0_chill <- round(par[2])
  T_base <- par[3]
  C_a <- par[4]
  C_b <- par[5]
  C_c <- par[6]
  C_ini <- par[7]
  F_crit <- par[8]
  C_req <- par[9]

  # chilling
  Rc <- 1 / ( 1 + exp( C_a * (data$Ti - C_c) ^ 2 + C_b * (data$Ti - C_c) ))
  Rc[1:t0_chill,] <- 0
  Sc <- apply(Rc,2, cumsum)

  k <- C_ini + Sc * (1 - C_ini)/C_req
  k[Sc >= C_req] <- 1

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Parallel M1 model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- PM1(data = data, par = par)
#'}

PM1 <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 8){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 <- par[1]
  T_base <- par[2]
  T_opt <- par[3]
  T_min <- par[4]
  T_max <- par[5]
  C_ini <- par[6]
  F_crit <- par[7]
  C_req <- par[8]

  # chilling
  Rc <- triangular_temperature_response(data$Ti,
                                       T_opt = T_opt,
                                       T_min = T_min,
                                       T_max = T_max)
  Rc[1:t0,] <- 0
  Rc[t0:nrow(Rc),] <- 0
  Sc <- apply(Rc,2, cumsum)

  k <- C_ini + Sc * (1 - C_ini)/C_req
  k[Sc >= C_req] <- 1

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- ((data$Li / 10) ^ k) * Rf
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Parallel M1 model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' using a bell shaped chilling response
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- PM1(data = data, par = par)
#'}

PM1b <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 8){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 = round(par[1])
  T_base = par[2]
  C_a = par[3]
  C_b = par[4]
  C_c = par[5]
  C_ini = par[6]
  F_crit = par[7]
  C_req = par[8]

  # chilling
  Rc = 1 / ( 1 + exp( C_a * (data$Ti - C_c) ^ 2 + C_b * (data$Ti - C_c) ))
  Rc[1:t0,] = 0
  Rc[t0:nrow(Rc),] = 0
  Sc = apply(Rc,2, cumsum)

  k = C_ini + Sc * (1 - C_ini)/C_req
  k[Sc >= C_req] = 1

  # forcing
  Rf = data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf = ((data$Li / 10) ^ k) * Rf
  Rf[1:t0,] = 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' PhotoThermal Time model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- PTT(data = data, par = par)
#'}

PTT <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1]) # int
  T_base <- par[2]
  F_crit <- par[3]

  # create forcing/chilling rate vector
  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- (data$Li / 24) * Rf
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' PhotoThermal Time model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' with a sigmoidal temperature response (Kramer 1994)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- PTTs(data = data, par = par)
#'}

PTTs <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- par[1]
  b <- par[2]
  c <- par[3]
  F_crit <- par[4]

  # create forcing/chilling rate vector
  # forcing
  Rf <- 1 / (1 + exp(-b * (data$Ti - c)))
  Rf <- (data$Li / 24) * Rf
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Standard growing season index model
#'
#' as defined by Xin et al. 2015 (Rem. Sens. Env.)
#'
#' No clear accumulation start date was set in the above mentioned
#' manuscript, as such we assume a start date of 21th of Dec, or 1th of Jan.
#' depending on the offset used in the generated validation data.
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- AGSI(data = data, par = par)
#'}

SGSI <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # exit the routine if data is missing
  if (is.null(data$Tmini) |
      is.null(data$Tmaxi) |
      is.null(data$VPDi)  |
      is.null(data$Li)
  ){
    stop("Not all required driver data is available")
  }

  # set start of accumulation period
  t0 <- which(data$doy == -11)
  if (length(t0) == 0){
    t0 <- 1
  }

  # extract the parameter values from the
  # par argument for readability
  Tmmin <- par[1]
  Tmmax <- par[2]
  F_crit <- par[3]
  VPDmin <- par[4]
  VPDmax <- par[5]
  photo_min <- par[6]
  photo_max <- par[7]

  # rescaling all parameters between 0 - 1 for the
  # acceptable ranges, if outside these ranges
  # set to 1 (unity) or 0 respectively
  Tmin <- (data$Tmini - Tmmin)/(Tmmax - Tmmin)
  VPD <- (data$VPDi - VPDmin)/(VPDmax - VPDmin)
  photo <- (data$Li - photo_min)/(photo_max - photo_min)

  # set outliers to 1 or 0 (step function)
  Tmin[which(data$Tmini <= Tmmin)] <- 0
  Tmin[which(data$Tmini >= Tmmax)] <- 1

  VPD[which(data$VPDi >= VPDmax)] <- 0
  VPD[which(data$VPDi <= VPDmin)] <- 1

  photo[which(data$Li <= photo_min)] <- 0
  photo[which(data$Li >= photo_max)] <- 1

  # calculate the index for every value
  # but set values for time t0 to 0
  GSI <- Tmin * VPD * photo
  GSI[1:t0,] <- 0

  # DOY of budburst criterium as calculated
  # by cummulating the GSI hence AGSI
  doy <- apply(GSI,2, function(xt){
    data$doy[which(zoo::rollmean(xt, 21, na.pad = TRUE) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Sequential model (M1 variant)
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- SM1(data = data, par = par)
#'}

SM1 <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 8){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 <- par[1]
  t0_chill <- par[2]
  T_base <- par[3]
  T_opt <- par[4]
  T_min <- par[5]
  T_max <- par[6]
  F_crit <- par[7]
  C_req <- par[8]

  # sanity check
  if (t0 <= t0_chill){
    return(rep(NA,ncol(data$Ti)))
  }

  # chilling
  Rc <- triangular_temperature_response(data$Ti,
                                       T_opt = T_opt,
                                       T_min = T_min,
                                       T_max = T_max)
  Rc[1:t0_chill,] <- 0
  Rc[t0:nrow(Rc),] <- 0
  Sc <- apply(Rc,2, cumsum)

  # chilling requirement has to be met before
  # accumulation starts (binary choice)
  k <- as.numeric(Sc >= C_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- (data$Li / 24) ^ k * Rf
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Sequential model (M1 variant)
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' using a bell shaped chilling response
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- SM1b(data = data, par = par)
#'}

SM1b <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 8){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  t0_chill <- round(par[2])
  T_base <- par[3]
  C_a <- par[4]
  C_b <- par[5]
  C_c <- par[6]
  F_crit <- par[7]
  C_req <- par[8]

  # sanity check
  if (t0 <= t0_chill){
    return(rep(NA,ncol(data$Ti)))
  }

  # chilling
  Rc <- 1 / ( 1 + exp( C_a * (data$Ti - C_c) ^ 2 + C_b * (data$Ti - C_c) ))
  Rc[1:t0_chill,] <- 0
  Rc[t0:nrow(Rc),] <- 0
  Sc <- apply(Rc,2, cumsum)

  # chilling requirement has to be met before
  # accumulation starts (binary choice)
  k <- as.numeric(Sc >= C_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- ((data$Li / 24) ^ k) * Rf
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Sequential model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = SQ(data = data, par = par)
#'}

SQ <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 8){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  t0_chill <- round(par[2])
  T_base <- par[3]
  T_opt <- par[4]
  T_min <- par[5]
  T_max <- par[6]
  F_crit <- par[7]
  C_req <- par[8]

  # sanity check t0 always comes after t0_chill
  if (t0 <= t0_chill){
    return(
      shape_model_output(data = data,
                         doy = rep(NA, ncol(data$Ti)))
    )
  }

  # chilling
  Rc <- triangular_temperature_response(
    data$Ti,
    T_opt = T_opt,
    T_min = T_min,
    T_max = T_max)

  Rc[1:t0_chill,] <- 0
  Rc[t0:nrow(Rc),] <- 0
  Sc <- apply(Rc, 2, cumsum)

  # chilling requirement has to be met before
  # accumulation starts (binary choice)
  k <- as.numeric(Sc >= C_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0 # CHECK THIS IN LITERATURE

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Sequential model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' using a bell shaped chilling response
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- SQb(data = data, par = par)
#'}

SQb <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 8){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  t0_chill <- round(par[2])
  T_base <- par[3]
  C_a <- par[4]
  C_b <- par[5]
  C_c <- par[6]
  F_crit <- par[7]
  C_req <- par[8]

  # chilling
  Rc <- 1 / ( 1 + exp( C_a * (data$Ti - C_c) ^ 2 + C_b * (data$Ti - C_c) ))
  Rc[1:t0_chill,] <- 0
  Rc[t0:nrow(Rc),] <- 0
  Sc <- apply(Rc,2, cumsum)

  # chilling requirement has to be met before
  # accumulation starts (binary choice)
  k <- as.numeric(Sc >= C_req)

  # forcing (with bell shaped curve -- only for forcing not chilling?)
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Thermal Time model
#'
#' simple growing degree day model as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = TT(data = data, par = par)
#'}

TT <- function(par, data ){

  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]

  # simple degree day sum setup
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Thermal Time model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' with a sigmoidal temperature response (Kramer 1994)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- TTs(data = data, par = par)
#'}

TTs <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  b <- par[2]
  c <- par[3]
  F_crit <- par[4]

  # sigmoid temperature response
  Rf <- 1 / (1 + exp(-b * (data$Ti - c)))
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Unified M1 model
#'
#' as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- UM1(data = data, par = par)
#'}

UM1 <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 8){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 <- round(par[1])
  T_base <- par[2]
  T_opt <- par[3]
  T_min <- par[4]
  T_max <- par[5]
  f <- par[6]
  w <- par[7]
  C_req <- par[8]

  # chilling accumulation using the
  # bell shaped temperature response
  Rc <- triangular_temperature_response(data$Ti,
                                       T_opt = T_opt,
                                       T_min = T_min,
                                       T_max = T_max)
  Rc[1:t0,] <- 0
  Sc <- apply(Rc, 2, cumsum)

  # chilling requirement has to be met before
  # accumulation starts (binary choice)
  k <- as.numeric(Sc >= C_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- ((data$Li / 10) ^ k) * Rf
  Rf[1:t0,] <- 0
  Sf <- apply(Rf, 2, cumsum)

  # DOY meeting F_crit, subtract the forcing matrix
  # from the F_crit matrix in order to speed things up
  # only the location of transition from - to + is
  # claculated to estimate the transition dates
  Sfc <- Sf - (w * exp(f * Sc))

  doy <- apply(Sfc, 2, function(x){
    data$doy[which(x > 0)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Accumulated growing season index model
#'
#' as defined by Xin et al. 2015 (Rem. Sens. Env.)
#'
#' The starting point of accumulation is not clearly indicated
#' in the publication so we assume December 21.
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- AGSI(data = data, par = par)
#'}

AGSI <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop(sprintf("model parameter(s) out of range (too many, too few, %s provided)",
                 length(par)))
  }

  # exit the routine if data is missing
  if (is.null(data$Tmini) |
      is.null(data$Tmaxi) |
      is.null(data$VPDi)  |
      is.null(data$Li)
  ){
    stop("Not all required driver data is available")
  }

  # set start of accumulation period
  t0 <- which(data$doy == -11)
  if (length(t0)==0){
    t0 <- 1
  }

  # extract the parameter values from the
  # par argument for readability
  Tmmin <- par[1]
  Tmmax <- par[2]
  F_crit <- par[3]
  VPDmin <- par[4]
  VPDmax <- par[5]
  photo_min <- par[6]
  photo_max <- par[7]

  # rescaling all parameters between 0 - 1 for the
  # acceptable ranges, if outside these ranges
  # set to 1 (unity) or 0 respectively
  Tmin <- (data$Tmini - Tmmin)/(Tmmax - Tmmin)
  VPD <- (data$VPDi - VPDmin)/(VPDmax - VPDmin)
  photo <- (data$Li - photo_min)/(photo_max - photo_min)

  # set outliers to 1 or 0
  Tmin[which(data$Tmini <= Tmmin)] <- 0
  Tmin[which(data$Tmini >= Tmmax)] <- 1

  VPD[which(data$VPDi >= VPDmax)] <- 0
  VPD[which(data$VPDi <= VPDmin)] <- 1

  photo[which(data$Li <= photo_min)] <- 0
  photo[which(data$Li >= photo_max)] <- 1

  # calculate the index for every value
  # but set values for time t0 to 0
  GSI <- Tmin * VPD * photo
  GSI[1:t0,] <- 0

  # DOY of budburst criterium as calculated
  # by cummulating the GSI hence AGSI
  doy <- apply(GSI,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' CU chilling degree model
#'
#' Chen et al. 2017 (Journal of Ecology)
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model, tropical forest
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- CU(data = data, par = par)
#'}

CU <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  ni <- round(par[1])
  nd <- round(par[2])
  T_base <- par[3]
  F_crit <- par[4]

  # simple chilling degree day sum setup
  # with lagged response
  data$Ti[data$Ti > T_base] <- 0
  data$Ti <- abs(data$Ti)

  # calculate floral induction time series
  Fd <- apply(data$Ti, 2, function(x){
    Fi <- zoo::rollapply(x, ni, sum, align = "right", fill = 0, na.rm = TRUE)
    Fd <- c(rep(0,nd),Fi[1:(length(Fi) - nd)])
  })

  # DOY of budburst criterium
  doy <- apply(Fd, 2, function(xt){
    data$doy[which(xt >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

#' Unified model
#'
#' as defined in Chuine 2000
#'
#' @param data input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @param par a vector of parameter values, this is functions specific
#' @return raster or vector with estimated phenophase timing (in DOY)
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate <- UN(data = data, par = par)
#'}

UN <- function(par, data){

  # This is an effort to reproduce the Unified Model
  # of Chuine 2000 in full form (not simplified)

  # exit the routine if parameters are missing
  if (length(par) != 9){
    stop("model parameter(s) out of range (too many, too few)")
  }

  CF <- function(x, a_c, b_c, c_c){
    1/(1 + exp(a_c * (x - c_c)^2 + b_c * (x - c_c)) )
  }

  # extract the parameter values from the
  # par argument in a more human readable form
  tc <- round(par[1]) # doy until when to accumulate chilling days
  a_c <- par[2] # sigmoid function chilling parameter a
  b_c <- par[3] # sigmoid function chilling parameter b
  c_c <- par[4] # sigmoid function chilling parameter c
  b_f <- par[5] # sigmoid function forcing parameter b
  c_f <- par[6] # sigmoid function forcing parameter c
  w <- par[7] # F* parameter w
  k <- par[8] # F* parameter k
  C_req <- par[9] # Chilling degree threshold requirement
  # i.e. C* whatever it is called

  # chilling accumulation using the
  # triangular shaped temperature response
  # basically convert normal temperatures to
  # what is called "chilling units" in the paper
  Rc <- CF(x = data$Ti, a_c, b_c, c_c)

  # Set values  < 0 to 0 (shouldn't count)
  # set NA values to 0 (when the output of the
  # triangular response is "empty" i.e. NA)
  Rc[is.na(Rc)] <- 0
  Rc[Rc < 0] <- 0

  # accumulate the chilling units, this is a
  # cummulative sum so you have all values along
  # the time axis (this saves time / iterations)
  Sc <- apply(Rc, 2, cumsum)

  # chilling requirement has to be met before
  # accumulation starts (binary choice) basically
  # binary mask to be applied to the Forcing temperature
  # data (sets anything before C_req to 0)
  m <- apply(Sc >= C_req, 2, as.numeric)

  # calculates when (row number) C_req is met
  row_loc <- apply(m,2,function(x)which(x == 1)[1])

  # if any row_loc is NA (C_req not met) or all
  # of m == 0 (same deal / fallback to be sure)
  # skip the rest as you won't be able to set
  # C_tot which by default should be > C_req
  if(any(is.na(row_loc)) | all(m == 0)){
    return(
      shape_model_output(
        data = data,
        doy = rep(9999,ncol(Sc))
      )
    )
  }

  # if all columns have valid values, and the associated
  # time location (row value) check if the tc value which
  # determines the total chilling degree day accumulation
  # exceeds the maximum value, if not C_tot < C_req which
  # is not allowed the total is always equal to or greater
  # then C_req. Skip if the condition is not met
  if (tc < max(row_loc)) {
    return(
      shape_model_output(
        data = data,
        doy = rep(9999, ncol(Sc))
      )
    )
  }

  # if above conditions are met
  # select the row which defines C_tot
  # this is non-dynamic across all sites / years
  # (as far as I can deduce from the Chuine paper)
  C_tot <- Sc[tc,]

  # apply the unified CF function
  # with a parameter set to 0
  Rf <- CF(x = data$Ti, 0, b_f, c_f)

  # Apply the chilling mask to forcing
  # temperature values
  Rf <- Rf  * m
                   
  # cummulate the forcing values
  Sfc <- apply(Rf, 2, cumsum)

  # calculate the Forcing requirement
  # based upon the C_tot value
  F_req <- w * exp(k * C_tot)

  # Trap invalid F* values, no need
  # to waste additional cycles
  if(any(is.na(F_req)) | any(is.infinite(F_req))){
    return(
      shape_model_output(
        data = data,
        doy = rep(9999,ncol(Sc))
      )
    )
  }

  # take the difference between the
  # forcing matrix and one filled with
  # the required F* values, where it
  # exceeds 0 first is the day of
  # leaf development
  Sfc <- sweep(Sfc, 2, F_req, FUN="-")
  doy <- apply(Sfc, 2, function(x){
    data$doy[which(x > 0)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  return(
    shape_model_output(
      data = data,
      doy = doy)
  )
}

###############Precip models#####################

#### Water-time (WT) ####

WT <- function(par, data ){

  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  F_crit <- par[2]
  P_base<- par[3]

  # simple degree day sum set-up, but for precip
  Rw <- data$Pi - P_base
  Rw[Rw < 0] <- 0
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}


####Water-time sigmoidal (WTs)####

WTs <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  b <- par[2]
  c <- par[3]
  F_crit <- par[4]

  # sigmoid precip response
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Photo-water time (PWT)####

PWT <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1]) # int
  F_crit <- par[2]
  P_base <- par[3]

  # create precip rate vector
  # forcing
  Rw <- data$Pi - P_base
  Rw[Rw < 0] <- 0
  Rw <- (data$Li / 24) * Rw
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Photo-water time sigmoidal (PWTs)####

PWTs <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  b <- par[2]
  c <- par[3]
  F_crit <- par[4]

  # create precip rate vector
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw <- (data$Li / 24) * Rw
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####M1 with precip (M1W)####

M1W <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  k <- par[2]
  F_crit <- par[3]
  P_base <- par[4]

  # create precip rate vector
  Rw <- data$Pi - P_base
  Rw[Rw < 0] <- 0
  Rw <- ((data$Li / 10) ^ k) * Rw
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####M1 with precip sigmoidal (M1Ws)####

M1Ws <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  b <- par[2]
  c <- par[3]
  k <- par[4]
  F_crit <- par[5]

  # create precip rate vector
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw <- ((data$Li / 10) ^ k) * Rw
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Sequential (precip then temp) (SQW)####

SQW <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  P_base <- par[4]
  P_req <- par[5]


  # Precip requirement
  Rw <- data$Pi - P_base
  Rw[Rw < 0] <- 0
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  # precip requirement has to be met before
  # temp accumulation starts (binary choice)
  k <- as.numeric(Sw >= P_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Sequential reverse (temp then precip) (SQWr)####
#r for reverse

SQWr <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  P_base <- par[4]
  T_req <- par[5]

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf[1:t0,] <- 0

  Sf <- apply(Rf, 2, cumsum)

  # chilling requirement has to be met before
  # accumulation starts (binary choice)
  k <- as.numeric(Sf >= T_req)

  # Precip requirement
  Rw <- data$Pi - P_base
  Rw[Rw < 0] <- 0
  Rw <- Rw * k
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Parallel (precip and temp) (PAW)####

PAW <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  P_base <- par[4]
  P_ini <- par[5]
  P_req <- par[6]

  # Precip requirement
  Rw <- data$Pi - P_base
  Rw[Rw < 0] <- 0
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  ##Determine k
  k <- P_ini + Sw * (1 - P_ini)/P_req
  k[Sw >= P_req] = 1

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Sequential sigmoidal (precip then temp) (SQWs)####

SQWs <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  P_req <- par[6]

  # sigmoid precip response
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  # SM requirement has to be met before
  # temp accumulation starts (binary choice)
  k <- as.numeric(Sw >= P_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Sequential reverse sigmoidal (temp then precip) (SQWrs)####

SQWrs <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  T_req <- par[6]


  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf[1:t0,] <- 0

  Sf <- apply(Rf, 2, cumsum)

  # temp requirement has to be met before
  # SM accumulation starts (binary choice)
  k <- as.numeric(Sf >= T_req)

  # sigmoid SM response
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw <- Rw * k

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Parallel sigmoidal (precip and temp) (PAW)####

PAWs <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  P_ini <- par[6]
  P_req <- par[7]


  # sigmoid SM response
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  ##Determine k
  k <- P_ini + Sw * (1 - P_ini)/P_req
  k[Sw >= P_req] = 1

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Sequential no P-base (SQW_NoPbase)####

SQW_NoPbase <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  P_req <- par[4]


  # Precip requirement
  Rw <- data$Pi
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  # precip requirement has to be met before
  # temp accumulation starts (binary choice)
  k <- as.numeric(Sw >= P_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####SQW_cdd: precip resets with consecutive dry days (cdd)####

SQW_cdd <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  #parameters
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  P_req <- par[4]
  cdd_thres <- round(par[5])

  #Precip requirement
  Rw <- data$Pi
  Rw[1:t0,] <- 0

  ##Calculate cdd

  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > 0, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)

  #force days before t0 to be 0 so not counted in cdd
  k_cdd[1:t0,] <- 0

  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, stats::ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }

  #Pull out years
  col_names <- names(k_cdd)

  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output

  }

  ##Reset Sw when hit cdd_thres

  #change matrices to dataframes
  Rw_df <- data.frame(Rw)
  colnames(Rw_df) <- paste0(data$site, "_", data$year)

  cdd_matrix_df <- data.frame(cdd_matrix)
  colnames(cdd_matrix_df) <- paste0(data$site, "_", data$year)

  #function to reset Sw when hit cdd_thres
  cumsum_func <- function(col){
    as.matrix(stats::ave(Rw_df[[col]], cumsum(cdd_matrix_df[[col]] == cdd_thres), FUN = cumsum))
  }

  #Make empty matrix
  Sw <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- cumsum_func(col = year)
    Sw[,i] <- output

  }

  ## precip requirement has to be met before temp accumulation

  # assign output matrix
  k <- matrix(0, nrow = 365, ncol = length(col_names))
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)

  # assign value of 1 once precip requirement met
  for (i in 1:cols){
    for (j in 1:rows){
      if(Sw[j,i] >= P_req){
        k[j:rows,i] <- 1
        break
      }
    }
  }


  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}


####SQW_Tmin: temp resets when Tmin below threshold####

SQW_Tmin <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  #parameters
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  P_req <- par[4]
  T_thres <- par[5]


  #Precip requirement
  Rw <- data$Pi
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  # precip requirement has to be met before
  # temp accumulation starts (binary choice)
  k <- as.numeric(Sw >= P_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  #change matrices to dataframes
  Rf_df <- data.frame(Rf)
  colnames(Rf_df) <- paste0(data$site, "_", data$year)

  Tmin_df <- data.frame(data$Tmini)
  colnames(Tmin_df) <- paste0(data$site, "_", data$year)

  #function to reset Sf when hit T_thres
  T_thres_func <- function(col){
    as.matrix(stats::ave(Rf_df[[col]], cumsum(Tmin_df[[col]] <= T_thres), FUN = cumsum))
  }

  #Pull out years
  col_names <- names(Tmin_df)

  #Make empty matrix
  Sf <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- T_thres_func(col = year)
    Sf[,i] <- output

  }

  # DOY of budburst criterium
  doy <- apply(Sf, 2, function(xt){
    data$doy[which(xt >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)

}



####SQW_cdd_Tmin: combine cdd & Tmin models####

SQW_cdd_Tmin <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }

  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  P_req <- par[4]
  cdd_thres <- round(par[5])
  T_thres <- par[6]

  # Precip requirement
  Rw <- data$Pi
  Rw[1:t0,] <- 0

  ##Calculate cdd

  #If precip = 0 (a dry day), then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > 0, 0, 1))
  colnames(k_cdd)<- paste0(data$site, "_", data$year)

  #force days before t0 to be 0 so not counted in cdd
  k_cdd[1:t0,] <- 0

  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(
      k_cdd,
      stats::ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }

  #Pull out years
  col_names <- names(k_cdd)

  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output

  }


  ##Reset Sw when hit cdd_thres

  #change matrices to dataframes
  Rw_df <- data.frame(Rw)
  colnames(Rw_df) <- paste0(data$site, "_", data$year)

  cdd_matrix_df <- data.frame(cdd_matrix)
  colnames(cdd_matrix_df) <- paste0(data$site, "_", data$year)

  #function to reset Sw when hit cdd_thres
  cumsum_func <- function(col){
    as.matrix(stats::ave(Rw_df[[col]], cumsum(cdd_matrix_df[[col]] == cdd_thres), FUN = cumsum))
  }

  #Make empty matrix
  Sw <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- cumsum_func(col = year)
    Sw[,i] <- output

  }


  ##precip requirement has to be met before temp accumulation starts

  # assign output matrix
  k <- matrix(0, nrow = 365, ncol = length(col_names))
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)

  # assign value of 1 once precip requirement met
  for (i in 1:cols){
    for (j in 1:rows){
      if(Sw[j,i] >= P_req){
        k[j:rows,i] <- 1
        break
      }
    }
  }


  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  #change matrices to dataframes
  Rf_df <- data.frame(Rf)
  colnames(Rf_df) <- paste0(data$site, "_", data$year)

  Tmin_df <- data.frame(data$Tmini)
  colnames(Tmin_df) <- paste0(data$site, "_", data$year)

  #function to reset Sf when hit T_thres
  T_thres_func <- function(col){
    as.matrix(stats::ave(Rf_df[[col]], cumsum(Tmin_df[[col]] <= T_thres), FUN = cumsum))
  }

  #Make empty matrix
  Sf <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- T_thres_func(col = year)
    Sf[,i] <- output

  }

  # DOY of budburst criterium
  doy <- apply(Sf, 2, function(xt){
    data$doy[which(xt >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)

}



####SQW_Pi_Tmin: precip reset with Tmin below threshold####

SQW_Pi_Tmin <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  #parameters
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  P_req <- par[4]
  T_thres <- par[5]

  # Precip requirement
  Rw <- data$Pi
  Rw[1:t0,] <- 0

  #change matrices to dataframes
  Rw_df <- data.frame(Rw)
  colnames(Rw_df) <- paste0(data$site, "_", data$year)

  Tmin_df <- data.frame(data$Tmini)
  colnames(Tmin_df) <- paste0(data$site, "_", data$year)

  #function to reset Sw when hit T_thres
  T_thres_func_Sw <- function(col){
    as.matrix(stats::ave(Rw_df[[col]], cumsum(Tmin_df[[col]] <= T_thres), FUN = cumsum))
  }

  #Pull out years
  col_names <- names(Tmin_df)

  #Make empty matrix
  Sw <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- T_thres_func_Sw(col = year)
    Sw[,i] <- output

  }

  ## precip requirement has to be met before temp accumulation starts

  # assign output matrix
  k <- matrix(0, nrow = 365, ncol = length(col_names))
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)

  # assign value of 1 once precip requirement met
  for (i in 1:cols){
    for (j in 1:rows){
      if(Sw[j,i] >= P_req){
        k[j:rows,i] <- 1
        break
      }
    }
  }

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)

}



####SQWs_cdd####

SQWs_cdd <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop("model parameter(s) out of range (too many, too few)")
  }

  #Parameters
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  P_req <- par[6]
  cdd_thres <- round(par[7])

  #Precip matrix to calculte cdd
  Precip <- data$Pi
  Precip[1:t0,] <- 0

  ##Calculate cdd

  #If precip = 0(a dry day), then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Precip > 0, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)

  #force days before t0 to be 0 so not counted in cdd
  k_cdd[1:t0,] <- 0

  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, stats::ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }

  #Pull out years
  col_names <- names(k_cdd)

  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output

  }


  #Calculate Rw - sigmoidal precip response
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw[1:t0,] <- 0


  ##Reset Sw when hit cdd_thres

  #change matrices to dataframes
  Rw_df <- data.frame(Rw)
  colnames(Rw_df) <-  paste0(data$site, "_", data$year)

  cdd_matrix_df <- data.frame(cdd_matrix)
  colnames(cdd_matrix_df) <-  paste0(data$site, "_", data$year)

  #function to reset Sw when hit cdd_thres
  cumsum_func <- function(col){
    as.matrix(stats::ave(Rw_df[[col]], cumsum(cdd_matrix_df[[col]] == cdd_thres), FUN = cumsum))
  }

  #Make empty matrix
  Sw <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- cumsum_func(col = year)
    Sw[,i] <- output

  }

  # precip requirement has to be met before temp accumulation starts

  # assign output matrix
  k <- matrix(0, nrow = 365, ncol = length(col_names))
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)

  # assign value of 1 once precip requirement met
  for (i in 1:cols){
    for (j in 1:rows){
      if(Sw[j,i] >= P_req){
        k[j:rows,i] <- 1
        break
      }
    }
  }

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer or a vector
  shape_model_output(data = data, doy = doy)
}


####SQWs_Tmin####

SQWs_Tmin <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop("model parameter(s) out of range (too many, too few)")
  }

  #Parameters
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  P_req <- par[6]
  T_thres <- par[7]

  #Calculate Rw - sigmoid precip response
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  # precip requirement has to be met before
  # temp accumulation starts (binary choice)
  k <- as.numeric(Sw >= P_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  #change matrices to dataframes
  Rf_df <- data.frame(Rf)
  colnames(Rf_df) <- paste0(data$site, "_", data$year)

  Tmin_df <- data.frame(data$Tmini)
  colnames(Tmin_df) <- paste0(data$site, "_", data$year)


  #function to reset Sf when hit T_thres
  T_thres_func <- function(col){
    as.matrix(stats::ave(Rf_df[[col]], cumsum(Tmin_df[[col]] <= T_thres), FUN = cumsum))
  }

  #Pull out years
  col_names <- names(Tmin_df)

  #Make empty matrix
  Sf <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- T_thres_func(col = year)
    Sf[,i] <- output

  }


  # DOY of budburst criterium
  doy <- apply(Sf, 2, function(xt){
    data$doy[which(xt >= F_crit)[1]]
  })

  # set export format, either a rasterLayer or a vector
  shape_model_output(data = data, doy = doy)

}



####SQWs_cdd_Tmin####

SQWs_cdd_Tmin <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 8){
    stop("model parameter(s) out of range (too many, too few)")
  }

  #Parameters
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  P_req <- par[6]
  cdd_thres <- round(par[7])
  T_thres <- par[8]

  #Precip matrix to calculte cdd
  Precip <- data$Pi
  Precip[1:t0,] <- 0

  ##Calculate cdd

  #If precip = 0(a dry day), then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Precip > 0, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)

  #force days before t0 to be 0 so not counted in cdd
  k_cdd[1:t0,] <- 0

  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, stats::ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }

  #Pull out years
  col_names <- names(k_cdd)

  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output

  }

  #Calculate Rw - sigmoidal precip response
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw[1:t0,] <- 0

  ##Reset Sw when hit cdd_thres

  #change matrices to dataframes
  Rw_df <- data.frame(Rw)
  colnames(Rw_df) <-  paste0(data$site, "_", data$year)

  cdd_matrix_df <- data.frame(cdd_matrix)
  colnames(cdd_matrix_df) <-  paste0(data$site, "_", data$year)

  #function to reset Sw when hit cdd_thres
  cumsum_func <- function(col){
    as.matrix(stats::ave(Rw_df[[col]], cumsum(cdd_matrix_df[[col]] == cdd_thres), FUN = cumsum))
  }

  #Make empty matrix
  Sw <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- cumsum_func(col = year)
    Sw[,i] <- output

  }

  # precip requirement has to be met before temp accumulation starts

  # assign output matrix
  k <- matrix(0, nrow = 365, ncol = length(col_names))
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)

  # assign value of 1 once precip requirement met
  for (i in 1:cols){
    for (j in 1:rows){
      if(Sw[j,i] >= P_req){
        k[j:rows,i] <- 1
        break
      }
    }
  }

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  #change matrices to dataframes
  Rf_df <- data.frame(Rf)
  colnames(Rf_df) <- paste0(data$site, "_", data$year)

  Tmin_df <- data.frame(data$Tmini)
  colnames(Tmin_df) <- paste0(data$site, "_", data$year)


  #function to reset Sf when hit T_thres
  T_thres_func <- function(col){
    as.matrix(stats::ave(Rf_df[[col]], cumsum(Tmin_df[[col]] <= T_thres), FUN = cumsum))
  }

  #Pull out years
  col_names <- names(Tmin_df)

  #Make empty matrix
  Sf <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- T_thres_func(col = year)
    Sf[,i] <- output

  }

  # DOY of budburst criterium
  doy <- apply(Sf, 2, function(xt){
    data$doy[which(xt >= F_crit)[1]]
  })

  # set export format, either a rasterLayer or a vector
  shape_model_output(data = data, doy = doy)

}



####SQWs_Pi_Tmin####

SQWs_Pi_Tmin <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop("model parameter(s) out of range (too many, too few)")
  }

  #parameters
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  P_req <- par[6]
  T_thres <- par[7]

  #Calculate Rw - sigmoid precip response
  Rw <- 1 / (1 + exp(-b * (data$Pi - c)))
  Rw[1:t0,] <- 0

  #change matrices to dataframes
  Rw_df <- data.frame(Rw)
  colnames(Rw_df) <- paste0(data$site, "_", data$year)

  Tmin_df <- data.frame(data$Tmini)
  colnames(Tmin_df) <- paste0(data$site, "_", data$year)

  #function to reset Sw when hit T_thres
  T_thres_func_Sw <- function(col){
    as.matrix(stats::ave(Rw_df[[col]], cumsum(Tmin_df[[col]] <= T_thres), FUN = cumsum))
  }

  #Pull out years
  col_names <- names(Tmin_df)

  #Make empty matrix
  Sw <- matrix(NA, nrow = 365, ncol = length(col_names))

  #loop years (columns) through function
  for (i in 1:length(col_names)) {

    year <- col_names[i]
    output <- T_thres_func_Sw(col = year)
    Sw[,i] <- output

  }

  ## precip requirement has to be met before temp accumulation starts

  # assign output matrix
  k <- matrix(0, nrow = 365, ncol = length(col_names))
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)

  # assign value of 1 once precip requirement met
  for (i in 1:cols){
    for (j in 1:rows){
      if(Sw[j,i] >= P_req){
        k[j:rows,i] <- 1
        break
      }
    }
  }

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)

}


############Soil moisture models##################


####Water Time (WT_SM)####

WT_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  F_crit <- par[2]
  SM_base <- par[3]

  # simple degree day sum set-up, but for SM
  Rw <- data$SM - SM_base
  Rw[Rw < 0] <- 0
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Photo-water time (PWT_SM)####

PWT_SM <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1]) # int
  F_crit <- par[2]
  SM_base <- par[3]

  # create SM rate vector
  # forcing
  Rw <- data$SM - SM_base
  Rw[Rw < 0] <- 0
  Rw <- (data$Li / 24) * Rw
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####M1 with SM (M1W_SM)####

M1W_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- par[1]
  k <- par[2]
  F_crit <- par[3]
  SM_base <- par[4]

  # create SM rate vector
  Rw <- data$SM - SM_base
  Rw[Rw < 0] <- 0
  Rw <- ((data$Li / 10) ^ k) * Rw
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Sequential (SM then temp) (SQW_SM)####

SQW_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  SM_base <- par[4]
  SM_req <- par[5]

  # SM requirement
  Rw <- data$SM - SM_base
  Rw[Rw < 0] <- 0
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  # SM requirement has to be met before
  # temp accumulation starts (binary choice)
  k <- as.numeric(Sw >= SM_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}


####Sequential reverse (temp then SM) (SQWr_SM)####

SQWr_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  T_req <- par[4]
  SM_base <- par[5]

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf[1:t0,] <- 0

  Sf <- apply(Rf, 2, cumsum)

  # chilling requirement has to be met before
  # accumulation starts (binary choice)
  k <- as.numeric(Sf >= T_req)

  # SM requirement
  Rw <- data$SM - SM_base
  Rw[Rw < 0] <- 0
  Rw <- Rw * k
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Parallel (SM and temp) (PAW_SM)####

PAW_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 <- round(par[1])
  T_base <- par[2]
  F_crit <- par[3]
  SM_base <- par[4]
  SM_req <- par[5]
  SM_ini <- par[6]

  # SM requirement
  Rw <- data$SM - SM_base
  Rw[Rw < 0] <- 0
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  ##Determine k
  k <- SM_ini + Sw * (1 - SM_ini)/SM_req
  k[Sw >= SM_req] = 1

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Water time sigmoidal (WTs_SM)####

WTs_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  b <- par[2]
  c <- par[3]
  F_crit <- par[4]

  # sigmoid SM response
  Rw <- 1 / (1 + exp(-b * (data$SM - c)))
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}


####Photo-water time sigmoidal (PWTs_SM)####

PWTs_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- par[1]
  b <- par[2]
  c <- par[3]
  F_crit <- par[4]

  # create SM rate vector
  Rw <- 1 / (1 + exp(-b * (data$SM - c)))
  Rw <- (data$Li / 24) * Rw
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####M1 with SM sigmoidal (M1Ws_SM)####

M1Ws_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- par[1]
  b <- par[2]
  c <- par[3]
  k <- par[4]
  F_crit <- par[5]

  # create SM rate vector
  Rw <- 1 / (1 + exp(-b * (data$SM - c)))
  Rw <- ((data$Li / 10) ^ k) * Rw
  Rw[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Sequential sigmoidal (SM then temp) (SQWs_SM)####

SQWs_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  SM_req <- par[6]

  # sigmoid SM response
  Rw <- 1 / (1 + exp(-b * (data$SM - c)))
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  # SM requirement has to be met before
  # temp accumulation starts (binary choice)
  k <- as.numeric(Sw >= SM_req)

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Sequential reverse sigmoidal (temp then SM) (SQWrs_SM)####

SQWrs_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  T_req <- par[6]


  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf[1:t0,] <- 0

  Sf <- apply(Rf, 2, cumsum)

  # temp requirement has to be met before
  # SM accumulation starts (binary choice)
  k <- as.numeric(Sf >= T_req)

  # sigmoid SM response
  Rw <- 1 / (1 + exp(-b * (data$SM - c)))
  Rw <- Rw * k

  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####Parallel sigmoidal (SM and temp) (PAWs_SM)####

PAWs_SM <- function(par, data){

  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument
  t0 <- round(par[1])
  T_base <- par[2]
  b <- par[3]
  c <- par[4]
  F_crit <- par[5]
  SM_req <- par[6]
  SM_ini <- par[7]


  # sigmoid SM response
  Rw <- 1 / (1 + exp(-b * (data$SM - c)))
  Rw[1:t0,] <- 0

  Sw <- apply(Rw, 2, cumsum)

  ##Determine k
  k <- SM_ini + Sw * (1 - SM_ini)/SM_req
  k[Sw >= SM_req] = 1

  # forcing
  Rf <- data$Ti - T_base
  Rf[Rf < 0] <- 0
  Rf <- Rf * k
  Rf[1:t0,] <- 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })

  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}


###############Autumn models#####################
#From Post & Richardson, 2025

####CDDP####
#Chilling Degree Day model with photoperiod adapted from PTT
#Existing model from Schadel et al., 2023

CDDP <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1]) # int
  T_base <- par[2]
  F_crit <- par[3]
  
  # create forcing/chilling rate vector
  # forcing
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  Rf <- ((1 - data$Li) / 24) * Rf
  Rf[1:t0,] <- 0
  
  # DOY of fall senescence criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}


####CDDs####
#Chilling Degree Day (CDD) adapted from
#with a sigmoidal temperature response
#Existing model from Schadel et al., 2023

CDDs <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  b <- par[2]
  c <- par[3]
  F_crit <- par[4]
  
  # sigmoid temperature response
  Rf <- 1 / (1 + exp(-b * (data$Tmini - c)))
  Rf[1:t0,] <- 0
  
  # DOY of fall senescence criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD2####
#t0 = SOS

CDD2 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 2){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base <- par[1]
  F_crit <- par[2]
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    Rf[1:SOS,i] <- 0
  }
  
  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####CDD3####
#t0 = SOS + dt

CDD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base <- par[1]
  F_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####CDDs3####
#t0 = SOS + dt

CDDs3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  b <- par[1]
  c <- par[2]
  F_crit <- par[3]
  dt <- round(par[4])
  
  # sigmoid temperature response
  Rf <- 1 / (1 + exp(-b * (data$Tmini - c)))
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  
  # DOY of fall senescence criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####SOS####
#EOS is certain number of days after SOS

SOS <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 1){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  dt <- round(par[1])
  
  #SOS plus # of days
  doy <- data$SOS_dates + dt
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####DD_W3####
#brown-down occurs after enough consecutive dry days (DD)
#t0 = SOS + dt

DD_W3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  P_crit <- par[1]
  P_base <- par[2]
  dt <- round(par[3])
  
  #Precip requirement
  Rw <- data$Pi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > P_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # DOY of budburst criterium
  doy <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= P_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####DD_SM3####
#t0 = SOS + dt

DD_SM3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  SM_base <- par[1]
  SM_crit <- par[2]
  dt <- round(par[3])
  
  #Precip requirement
  Rw <- data$SM
  
  ##Calculate cdd
  
  #If SM < SM_base, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > SM_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # DOY of budburst criterium
  doy <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= SM_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####DD_VPD3#### 
#t0 = SOS + dt

DD_VPD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  VPD_base <- par[1]
  VPD_crit <- par[2]
  dt <- round(par[3])
  
  #Precip requirement
  Rw <- data$VPDi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw < VPD_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when VPD above VPD base
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # DOY of budburst criterium
  doy <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= VPD_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}




####W3####
#t0 = SOS + dt

W3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  a <- par[1]
  P_crit <- par[2]
  dt <- round(par[3])
  dp <- round(par[4])
  
  #Make output matrix
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)
  Rw <- matrix(0,rows,cols)
  k <- c(1:floor(dp))
  W <- 0
  
  #Make precip matrix
  Pi_matrix <- data$Pi
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  # loop through to find precip over last "dp" days with decay function 
  
  for (i in 1:cols){
    
    for (j in 1:rows){
      
      t0 <- SOS_day[i] + dt
      
      if (j <= t0) {
        next
        
      }else{
        
        if(j > dp){
          
          (my_vector <- Pi_matrix[j-k,i]*(a^k))
          W <- sum(my_vector) + Pi_matrix[j,i]
        }
      }
      
      if (W <= P_crit){
        Rw[j,i] <- 1000
        break
      }
    }
  }
  
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(xt == 1000)[1]]
  })
  
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####SM3####
#t0 = SOS + dt

SM3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  SM_base <- par[1]
  SM_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rw <- data$SM - SM_base
  Rw[Rw > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) <= SM_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####VPD3####
#t0 = SOS + dt

VPD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  VPD_base <- par[1]
  VPD_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rw <- data$VPDi - VPD_base
  Rw[Rw < 0] <- 0
  
  #Set any VPD before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if(t0 <= 365){
      Rw[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= VPD_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####CDD_DD_W3####
#t0 = SOS + dt

CDD_DD_W3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  P_crit <- par[3]
  P_base <- par[4]
  dt <- round(par[5])
  
  #Precip requirement
  Rw <- data$Pi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > P_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= P_crit)[1]]})
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_DD_SM3####
#t0 = SOS + dt

CDD_DD_SM3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  SM_base <- par[3]
  SM_crit <- round(par[4])
  dt <- round(par[5])
  
  #Precip requirement
  Rw <- data$SM
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > SM_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= SM_crit)[1]]})
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_DD_VPD3####
#t0 = SOS + dt

CDD_DD_VPD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  VPD_base <- par[3]
  VPD_crit <- round(par[4])
  dt <- round(par[5])
  
  #Precip requirement
  Rw <- data$VPDi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw < VPD_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= VPD_crit)[1]]})
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_W3####
#t0 = SOS + dt

CDD_W3 <- function(par, data ){
  
  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  T_base <- par[1]
  a <- par[2]
  F_crit <- par[3]
  P_crit <- par[4]
  dt <- round(par[5])
  dp <- round(par[6])
  
  # create chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  
  #Calculate precip#
  
  #Make output matrix
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)
  Rw <- matrix(0,rows,cols)
  k <- c(1:floor(dp))
  W <- 0
  
  #Make precip matrix
  Pi_matrix <- data$Pi
  
  # loop through to find precip over last "dp" days with decay function 
  
  for (i in 1:cols){
    
    for (j in 1:rows){
      
      t0 <- SOS_day[i] + dt
      
      if (j <= t0) {
        next
        
      }else{
        
        if(j > dp){
          
          (my_vector <- Pi_matrix[j-k,i]*(a^k))
          W <- sum(my_vector) + Pi_matrix[j,i]
        }
      }
      
      if (W <= P_crit){
        Rw[j,i] <- 1000
        break
      }
    }
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(Rw,2, function(xt){
    data$doy[which(xt == 1000)[1]]
  })
  
  # DOY of budburst criterium
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_SM3####

CDD_SM3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  SM_base <- par[3]
  SM_crit <- par[4]
  dt <- round(par[5])
  
  # create VPD forcing rate vector
  Rw <- data$SM - SM_base
  Rw[Rw > 0] <- 0
  
  #Set any VPD before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  
  # create chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) <= SM_crit)[1]]
  })
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_VPD3####
#t0 = SOS + dt

CDD_VPD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  VPD_base <- par[3]
  VPD_crit <- par[4]
  dt <- round(par[5])
  
  # create VPD forcing rate vector
  Rw <- data$VPDi - VPD_base
  Rw[Rw < 0] <- 0
  
  #Set any VPD before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  
  # create chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= VPD_crit)[1]]
  })
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}


####W_T_bal3####
#t0 = SOS + dt


W_T_bal3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  a <- par[1]
  F_crit <- par[2]
  dt <- round(par[3])
  
  
  #Precip requirement
  Rw <- data$Pi - (a*data$Ti)
  
  #Set any precip before SOS to 20 (so not triggered)
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  
  # DOY of fall senescence criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDDP3####
#t0 = SOS + dt
#adjusted so same direction of influence on Rf for darker and colder

CDDP3 <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base <- par[1]
  F_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  # forcing
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  #temper Rf based on Li
  Rf <- -1 * ((1-(data$Li/24)) * Rf)
  
  # DOY of fall senescence criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####SMP3####
#t0 = SOS + dt

SMP3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  SM_base <- par[1]
  SM_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rw <- data$SM - SM_base
  Rw[Rw > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  #temper Rf based on Li
  Rw <- (1-(data$Li/24)) * Rw
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) <= SM_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####VPDP3####
#VPD & photoperiod

VPDP3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  VPD_base <- par[1]
  F_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rw <- data$VPD - VPD_base
  Rw[Rw < 0] <- 0
  
  #Set any VPD before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
    
  }
  
  
  #temper Rw based on Li
  Rw <- (1-(data$Li/24)) * Rw
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}


####CDDP_DD_W3####
#t0 = SOS + dt

CDDP_DD_W3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  P_crit <- par[3]
  P_base <- par[4]
  dt <- round(par[5])
  
  #Precip requirement
  Rw <- data$Pi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > P_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # create forcing/chilling rate vector
  # forcing
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  
  #Set any precip before SOS to 0
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  #temper Rf based on Li
  Rf <- -1 * ((1-(data$Li/24)) * Rf)
  
  # DOY of budburst criterium
  doy_1 <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= P_crit)[1]]})
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}

                   

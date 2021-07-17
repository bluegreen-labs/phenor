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
  Rc <- data$Ti - T_base
  Rc[Rc < 0] <- 1
  Rc[Rc >= 0] <- 0
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

  # Apply the chilling mask to forcing
  # temperature values, set anything < 0
  # to 0 (see )
  Rf <- data$Ti * m

  # apply the unified CF function
  # with a parameter set to 0
  Rf <- CF(x = Rf, 0, b_f, c_f)

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

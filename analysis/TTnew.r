#' Thermal Time model as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#'
#' @param data: a nested list of data
#' @param par: a vector of parameter values, this is functions specific
#' @keywords phenology, model
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = TT(data = data, par = par)
#'}

library(phenor)

TTnew = function(par, data, plot = FALSE){

  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop("model parameter(s) out of range (too many, too few)")
  }

  # extract the parameter values from the
  # par argument for readability
  t0 = round(par[1])
  T_base = par[2]
  a = par[3]
  b = par[4]
  c = par[5]
  F_crit = par[6]
  L_base = par[7]

  # t1
  t1 = which(data$Li[,1] == min(data$Li[,1]))
  if (t0 < t1){
    return(rep(9999,length(data$transition_dates)))
  }

  # calculate MAT (latitude dependent of sorts)
  T = apply(data$ltm, 2, mean)

  # photopheriod
  Li = data$Li
  Li = t(apply(Li, 1,
              function(x){
                (x / 24) / ((a/(1 + exp(b * (T - c) ))) + 1)
                }))

  # photopheriod
  # Li = data$Li
  # Li[1:t0,] = 0
  # t1 = apply(Li,2, function(x){
  #  ifelse(x >= L_base, 1, 0)
  # })

  #Rf = t(apply(data$Ti, 1,
  #            function(x){
  #             x - (T_base * 1/((a/(1 + exp(b * (T - c) ))) + 1))
  #              }))

  # simple degree day sum setup
  Rf = data$Ti - T_base
  Rf[Rf < 0] = 0
  Rf = Li * Rf
  Rf[1:t0,] = 0

  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = data$doy[which(cumsum(xt) >= F_crit)[1]]
    doy[is.na(doy)] = 9999
    return(doy)
  })

  return(doy)
}

par = model_validation(model = "TTnew",
                par_ranges = "/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/parameter_ranges.csv",
                control = list(max.call = 5000, temperature = 10000))

print(par)
df = flat_format(phenocam_DB)
predicted = TTnew(par = par, data = df, plot = TRUE)
measured = df$transition_dates
T = apply(df$ltm, 2, mean)
photoperiod = (par[3]/(1 + exp(par[4] * (T - par[5]) ))) + 1

index = 1:length(T)
pred = data.frame(measured,predicted, photoperiod,index, T)

p = ggplot(data = pred,
           aes(x = index,
               y = T,
               colour = photoperiod)) + geom_point() + theme_bw() + scale_colour_gradient(low = "red", high = "blue")
plot(p)

Li = df$Li
Li = t(apply(Li, 1,
             function(x){
               (x / 24) / ((par[3]/(1 + exp(par[4] * (T - par[5]) ))) + 1)
             }))

Ti = t(apply(df$Ti, 1,
             function(x){
               x - (par[2] * 1/((par[3]/(1 + exp(par[4] * (T - par[5]) ))) + 1))
             }))

plot(Li[,1], type = 'n', ylim = c(0,0.5))
colours = heat.colors(length(T))
for (i in 1:length(T)){
  loc = order(T)[i]
  lines(Li[,loc], col = colours[i])
}
plot(T[order(photoperiod)],sort(photoperiod))

plot(Ti[,1], type = 'n', ylim = c(0,100))
colours = heat.colors(length(T))
for (i in 1:length(T)){
  loc = order(T)[i]
  lines(Ti[,loc], col = colours[i])
}

p = ggplot(data = pred ,
           aes(x = measured,
               y = predicted,
               colour = T)) + geom_point() + theme_bw() +
  scale_colour_gradient2(low = "red", mid = "green", high = "blue", midpoint = mean(T))
plot(p)

function(par = par, data = data){

# assign parameters
b1      = par[1]
b2 		  = par[2]
b3 		  = par[3]
b4 		  = par[4]
L       = par[5]
Phmin   = par[6]
Topt		= par[8]
Phmax		= par[9]

# assign constant values to declared variables Vmin/Vmax
Vmin 	= 0.001	# set to a small value but not 0
Vmax 	= 1.		# set to 100% [GCC values are scaled between 0-1
d     = 0		 	# decay flag
                           
# assign values to both the W and V arrays
# these are initial conditions that influence the final
# results to a great extent##
W		     = 0.
Wstart   = 0.
V 	     = 0.001
Sd 	     = 0.
b1 	     = Wp
m	       = 3600

# set soil water capacity
Wcap = data$wcap

# temperature response function returns
# a matrix when provided one
g = function(T,Topt){
  Tmax 	= 45; Tmin  = 0
  g = ((Tmax - T)/(Tmax - Topt)) * (((T - Tmin) / (Topt - Tmin)) ^ (Topt/(Tmax - Topt)))
  g[is.nan(g)] = 0
  return(g)
}

# number of rows in the matrix 
values = nrow(data$Pi) - 1
sites = ncol(data$Pi)


for (j in 1:sites){
for (i in 1:values){ # -- main loop
                           
  # first if statement traps initial conditions
  # where the values would be out of bound [looking
  # at place such as -1]
  # within the statement define Dtl and Dtl1
  # state of the soil water at time day and day-1
  # including a lag of L days
  if ( i - L - 1 < 0 ){
    Dt[i] = max(0, W[i] - b1)
    Dtl = Wstart
    Dtl1 = Wstart
  } else {
    Dt[i] = max(0, W[i] - b1)
    Dtl = Dt[i-L]
    Dtl1 = Dt[i-L-1]
  }
                           
  # if there is more precipitation [SWC]
  # compared to the previous day then decay = 0
  # otherwise the decay parameter is set to 1
  # and senescence sets in
  if ( Dtl > Dtl1 ){
    d = 0
  } else {
    d = 1
  }
                           
  # set dormancy conditions based upon the
  # available radiation
  dor = (Ra[i] - Phmin)/(Phmax - Phmin)
  if ( Ra[i] >= Phmax ){
    dor = 1
  }
  if ( Ra[i] <= Phmin ){
    dor = 0
  }
                           
  # SOIL WATER CONTENT SECTION:
  W[i+1] = W[i] + precip[i] - (1-V[i]) * (( Dt[i]/(Wcap - b1))^2) * evap[i] - g * b4 * V[i] * Dt[i]
                           
  # don't allow negative SWC values
  W[i+1] = max(0.0,min(Wcap,W[i+1]))
                           
  # VEGETATION GROWTH SECTION:
  V[i+1] = V[i] + g * dor * b2 * Dtl * (1 - V[i]/Vmax) - d * b3 * V[i] * (1-V[i])
                           
  # do not allow negative vegetation cover or values > 1
  V[i+1] = max(Vmin,min(Vmax,V[i+1]))
}
}
}
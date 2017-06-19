# Estimation of difference between
# the beginning and the end of the century
# for a CMIP5 scenario.
#
# This code is stateless (with floating parameters)
# so slight differences will exists between reported
# and generated values.
#
# This code and the statistics are meant as a
# methodological example NOT a formal analysis.
# The latter we leave to users of the phenor package.

# load the necessary libraries
library(phenor)

# create a set of demo parameters
# optimize for default Deciduous Broadleaf
# data (included in the phenor package)
# I use the simplest Thermal Time model (TT)
par = model_validation(model = "TT",
                       control = list(max.call = 10000,
                                      temperature = 10000))


# download year 2011 and 2100 for the CMIP5 run
data_2011 = format_cmip5(year = 2011)
data_2100 = format_cmip5(year = 2100)

# calculate model output
map2011 = estimate_phenology(par = par,
                             data = data_2011)
map2100 = estimate_phenology(par = par,
                             data = data_2100)

# take the difference between the two maps
mapdiff = map2100 - map2011

# report the mean and sd for the onset dates
# as calculate for the IGBP 4 (deciduous broadleaf forest) class
cat(sprintf("%s +_ %s\n",
            mean(mapdiff[igbp_4 > 0.2]),
            sd(mapdiff[igbp_4 > 0.2])))

# report min and max difference values
cat(sprintf("max: %s; min: %s",
  max(mapdiff[igbp_4 > 0.2]),
  min(mapdiff[igbp_4 > 0.2])))

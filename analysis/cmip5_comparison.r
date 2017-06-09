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

library(phenor)

# quick estimate of parameters for the default TT model
par = model_validation()

# download year 2011 and 2100 for the CMIP5 run
data_2011 = format_cmip5(year = 2011)
data_2100 = format_cmip5(year = 2100)

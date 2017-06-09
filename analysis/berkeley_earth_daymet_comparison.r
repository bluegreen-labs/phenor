# Statistics on the effect of spatial aggregation
# to 1x1 degrees on phenology estimates. Here we
# compare 1x1 degree Berkeley Earth phenor model
# output to Daymet model output at 1km resolution
# (upscaled to 1x1 degree)

# load the necessary libraries
library(phenor)

# path to store phenor data
path = "~/phenor_files"

# set working directory, create path parameters in
# necessary functions (TODO)
setwd(path)

# create a set of demo parameters
# optimize for default Deciduous Broadleaf
# data (included in the phenor package)
# I use the simplest Thermal Time model (TT)
par = model_validation(model = "TT",
                       control = list(max.call = 10000,
                                      temperature = 10000))

# download berkeley earth data (global data but subset for CONUS)
# for the same year as the Daymet data (scale comparison)
be_data = format_berkeley_earth(year = 2011)

# create maps using estimated parameters and
# the estimate_phenology routine
be_map = estimate_phenology(par = par,
                            model = "TT",
                            data = be_data)

daymet_map = estimate_phenology(par = par,
                                model = "TT",
                                path = path)

# crop the Berkeley Earth map
be_map = crop(be_map, extent(daymet_map))

# upscale the daymet data
zones = be_map
zones[] = 1:prod(dim(be_map))

# resmaple to match the extent of the Daymet map
# use a nearest neighbour approach (retains original cell numbers)
zones = resample(zones, daymet_map, method = "ngb")

# aggregate by zones, output is not a map!
mean_daymet = zonal(daymet_map, zones, fun = "mean", na.rm = TRUE)

# combine the data from the 1x1 degree Berkeley map
# with the data from the downsampled
be_daymet_stats = na.omit(cbind(as.vector(as.vector(be_map)), mean_daymet[, 2]))

# print output of the MAE + SD
model_difference = apply(test,1,function(x)abs(x[1]-x[2]))
cat(sprintf("%s +_ %s",
            mean(model_difference),
            sd(model_difference)
            ))


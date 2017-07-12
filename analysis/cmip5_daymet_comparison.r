# Statistics on the effect of spatial aggregation
# to 1x1 degrees on phenology estimates. Here we
# compare 1x1 degree CMIP5 phenor model (1 model only no ensemble)
# output to Daymet model output at 1km resolution
# (upscaled to 1x1 degree)

# load the necessary libraries
library(phenor)

# path to store phenor data
path = "~/phenor_files"

# set working directory, create path parameters in
# necessary functions (TODO)
setwd(path)

# download daymet tiles and format the data
# this might take a while (be warned), first
# list all the tiles needed
tiles = c(
  11936, 12295, 12296, 12297,
  12114, 12115, 12116, 12117,
  11934, 11935, 11754, 11755,
  11756, 12294, 11937)

# download all daymet tiles, include the preceeding year
# as we need this data too (see start_yr).
daymetr::download_daymet_tiles(tiles = tiles,
                               param = c("tmin","tmax"),
                               start = 2010,
                               end = 2011)

# now calculate the mean daily temperature from tmin and tmax
lapply(tiles, function(x)daymetr::daymet_tmean(path = path,
                                               tile = x,
                                               year = 2010,
                                               internal = FALSE))

lapply(tiles, function(x)daymetr::daymet_tmean(path = path,
                                               tile = x,
                                               year = 2011,
                                               internal = FALSE))

# now format all tiles acoording to the phenor format
# (default settings)
lapply(tiles, function(x)daymet_data = format_daymet_tiles(
  path = path,
  tile = x,
  year = 2011,
  internal = FALSE))

# create a set of demo parameters
# optimize for default Deciduous Broadleaf
# data (included in the phenor package)
# I use the simplest Thermal Time model (TT)
par = model_validation(model = "TT",
                       control = list(max.call = 10000,
                                      temperature = 10000))

# download CMIP5 data for the same year as the Daymet data
data_2011 = format_cmip5(year = 2011)

# create maps using estimated parameters and
# the estimate_phenology routine
cmip5_map_2011 = estimate_phenology(par = par,
                            model = "TT",
                            data = cmip5_data_2011)

daymet_map = estimate_phenology(par = par,
                                model = "TT",
                                path = path)

# crop the Berkeley Earth map
cmip5_map_2011 = crop(cmip5_map_2011, extent(daymet_map))

# upscale the daymet data
zones = cmip5_map_2011
zones[] = 1:prod(dim(cmip5_map_2011))

# resmaple to match the extent of the Daymet map
# use a nearest neighbour approach (retains original cell numbers)
zones = resample(zones, daymet_map, method = "ngb")

# aggregate by zones, output is not a map!
mean_daymet = zonal(daymet_map, zones, fun = "mean", na.rm = TRUE)

# combine the data from the 1x1 degree Berkeley map
# with the data from the downsampled
daymet_stats = na.omit(cbind(as.vector(as.vector(cmip5_map_2011)), mean_daymet[, 2]))

# print output of the MAE + SD
model_difference = apply(daymet_stats,1,function(x)abs(diff(x)))
cat(sprintf("%s +_ %s",
            mean(model_difference),
            sd(model_difference)
            ))


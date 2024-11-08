library(ecmwfr)
library(daymetr)
library(terra)

source("R/pr_fm_daymet_tiles.R")
source("R/pr_flatten.R")
source("R/pr_predict.R")
source("R/phenology_models.R")
source("R/helper_functions.R")

# download_daymet_tiles(
#   tiles = 11935,
#   start = 1980,
#   end = 1981,
#   param = c("tmin","tmax"),
#   path = "~/Downloads/"
# )

# daymet_grid_tmean(
#   path = "~/Downloads/",
#   product = 11935,
#   year = 1980,
#   internal = FALSE
# )
#
# daymet_grid_tmean(
#   path = "~/Downloads/",
#   product = 11935,
#   year = 1981,
#   internal = FALSE
# )

# drivers <- pr_fm_daymet_tiles(
#   path = "~/Downloads/",
#   tile = 11935,
#   year = 1981,
#   offset = 264,
#   internal = TRUE
# )
#
# set.seed(1234)
# optim.par <- pr_fit(
#   data = phenocam_DB,
#   cost = rmse,
#   model = "TT",
#   method = "GenSA"
# )

output <- pr_predict(
  optim.par$par,
  data = drivers,
  model = "TT"
  )

terra::plot(output)

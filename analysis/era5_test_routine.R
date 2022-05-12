#library(phenor)
# library(ecmwfr)
# options(keyring_backend="file")
# source("R/pr_dl_era5.R")
#
# data <- pr_dl_era5(
#   path = "~/Desktop",
#   user = "2088",
#   extent = c(
#     45,
#     -4,
#     42,
#     1.2
#     )
#   )
#
# r <- terra::rast("~/Desktop/era5.nc")
# terra::plot(r$t2m_1)
# maps::map("world", add = TRUE)
library(phenor)
source("R/pr_fm_era5.R")

data <- pr_fm_era5(
  path = "~/Desktop/",
  year = 2019,
  extent = c(
      -4,
      1.1,
      41,
      45
      )
)


# load the included data using
data("phenocam_DB")

# optimize model parameters
set.seed(1234)
optim.par <- pr_fit(
  data = phenocam_DB,
  cost = rmse,
  model = "TT",
  method = "GenSA"
)

output <- pr_predict(
  optim.par$par,
  data = data,
  model = "TT"
  )

raster::plot(output)

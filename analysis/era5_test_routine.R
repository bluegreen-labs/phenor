library(phenor)
library(ecmwfr)
options(keyring_backend="file")

source("R/pr_dl_era5.R")

# ERA5-land example
pr_dl_era5(
  path = "~/Desktop",
  user = "2088",
  product = "land",
  file = "era5-land.nc",
  extent = c(
    48.345183009475015,
    4.782986565238425,
    45.64153500799514,
    11.44515031394129
  )
)

source("R/pr_fm_era5.R")

data <- pr_fm_era5(
  path = "~/Desktop/",
  file = "era5-land.nc",
  year = 2019
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
maps::map("world", add = TRUE)

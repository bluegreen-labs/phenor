library(phenor)
library(ecmwfr)
options(keyring_backend="file")

source("R/pr_dl_era5.R")

data <- pr_dl_era5(
  path = "~/Desktop",
  user = "2088",
  product = "land",
  extent = c(
    -44.73496340744479,
    -78.4726872352846,
    -55.766660401835935,
    -58.49660835383497
    )
  )

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

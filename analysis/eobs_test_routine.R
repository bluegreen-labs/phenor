library(phenor)

source("R/pr_fm_eobs.R")

data <- pr_fm_eobs(
  path = "~/Desktop/",
  file = "xxxx",
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

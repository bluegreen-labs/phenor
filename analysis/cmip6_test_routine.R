library(phenor)
source("R/pr_fm_cmip.R")
source("R/pr_dl_cmip.R")

# format the data (year 2050)
data <- pr_dl_cmip(
  path = "~/Desktop",
  user = "2088",
  extent = c(
    50.73149477111302,
    -7.08887567473501,
    40.365567456020266,
    12.748594284373073
  )
)

data <- try(
  pr_fm_cmip(
    year = 2100,
    internal = TRUE,
    extent = c(
          50.73149477111302,
          -7.08887567473501,
          40.365567456020266,
          12.748594284373073
    ),
    path = "~/Desktop/"
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
maps::map("world", add = TRUE)

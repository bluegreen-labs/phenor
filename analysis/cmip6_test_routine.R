#library(phenor)

path <- system.file("extdata", package = "phenor")

# format the data (year 2050)
source("R/pr_fm_cmip.R")
data <- try(pr_fm_cmip(year = 2050, internal = TRUE, path = root))


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

bla <- pr_predict(optim.par$par, data = data, model = "TT")

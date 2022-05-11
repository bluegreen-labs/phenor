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

source("R/pr_fm_era5.R")

pr_fm_era5(
  path = "~/Desktop/",
  year = 2019,
    extent = c(
      45,
      -4,
      42,
      1.2
      )
)

#library(phenor)
#options(keyring_backend="file")

# download data (MIROC6 model 2000 - 2100)
#pr_dl_cmip(user = "2088")

# format the data (year 2050)
source("R/pr_fm_cmip.R")
data <- try(pr_fm_cmip(year = 2050, internal = TRUE))

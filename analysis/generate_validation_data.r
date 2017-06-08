# Generates validation data sets

# pull data from github
#system("git clone phenocamdata")

# move into directory
#setwd("~/phenocam_data/data_set_6")

# save data in the phenor package location
setwd("/data/Dropbox/Research_Projects/code_repository/bitbucket/PhenoR/data/")

# generate DB dataset
phenocam_DB = format_phenocam(path = "/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/data_record_6/DB",
                              threshold = 25)
devtools::use_data(phenocam_DB, overwrite = TRUE)

# generate EN dataset
phenocam_EN = format_phenocam(path = "/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/data_record_6/EN",
                              threshold = 25)
devtools::use_data(phenocam_EN, overwrite = TRUE)

# generate GR dataset (removed kendall, ibp, coaloilpoint)
phenocam_GR = format_phenocam(path = "/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/data_record_6/GR",
                              threshold = 25)
devtools::use_data(phenocam_GR, overwrite = TRUE)

# commit igbp masks to the repo
igbp_1 = raster::raster("../inst/extdata/igbp_1.tif")
devtools::use_data(igbp_1, overwrite = TRUE)

# commit igbp masks to the repo
igbp_4 = raster::raster("../inst/extdata/igbp_4.tif")
devtools::use_data(igbp_4, overwrite = TRUE)

# commit igbp masks to the repo
igbp_5 = raster::raster("../inst/extdata/igbp_5.tif")
devtools::use_data(igbp_5, overwrite = TRUE)

# commit igbp masks to the repo
igbp_10 = raster::raster("../inst/extdata/igbp_10.tif")
devtools::use_data(igbp_10, overwrite = TRUE)

# Generates validation data sets
# this for internal development use only.

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

# commit igbp masks to the repo / USE brick() otherwise only a
# reference to the file is written to file. With brick() data stored in
# memory is written to file.
# igbp_1 = raster::brick("/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/modis_data/igbp_1.tif")
# devtools::use_data(igbp_1, overwrite = TRUE)
#
# # commit igbp masks to the repo
# igbp_4 = raster::brick("/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/modis_data/igbp_4.tif")
# devtools::use_data(igbp_4, overwrite = TRUE)
#
# # commit igbp masks to the repo
# igbp_5 = raster::brick("/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/modis_data/igbp_5.tif")
# devtools::use_data(igbp_5, overwrite = TRUE)
#
# # commit igbp masks to the repo
# igbp_10 = raster::brick("/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/modis_data/igbp_10.tif")
# devtools::use_data(igbp_10, overwrite = TRUE)

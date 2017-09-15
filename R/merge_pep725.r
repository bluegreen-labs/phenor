#' Preprocessing of PEP725 data, merges separate files into tidy
#' data, with each observation a line, each column a different
#' parameter value.
#'
#' @param path: a path to the PEP725 data (species files only)
#' @return concatted data of all data in the path as a tidy data frame
#' listing PEP_ID, BBCH, YEAR, DAY, National_ID, LON, LAT, ALT, NAME
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' tidy_pep_data = merge_pep725()
#'}

merge_pep725 = function(path = "~"){

  # list all files
  pep_files = list.files(path, "*.csv", full.names = TRUE)

  # list file type locations
  data_file_locations = !grepl("^.*PEP725_BBCH.csv$",pep_files)
  station_file_locations = grepl("^.*stations\\.csv$",pep_files)

  # extract only the true data files and station info files
  # drop the BBCH and README data (but don't delete it - delist)
  data_files = pep_files[data_file_locations & !station_file_locations]
  station_files = pep_files[station_file_locations]

  # merge all the csv files for both the true observations and the
  # station meta-data (name, lat, lon etc.)
  observation_data = do.call(rbind,
                 lapply(data_files, function(file){
                   utils::read.csv2(file, sep = ";")
                   })
                 )

  station_locations = do.call(rbind,
                     lapply(station_files, function(file){
                      tmp = utils::read.csv2(file,
                                        sep = ";",
                                        stringsAsFactors = FALSE)
                      tmp$LON = as.numeric(tmp$LON)
                      tmp$LAT = as.numeric(tmp$LAT)
                      return(tmp)
                     })
  )

  # do a left merge to combine the observational data and the
  # station location meta-data returning basically the original
  # database structure (as a tidy file)
  pep_data = merge(observation_data, station_locations, by = "PEP_ID")

  # return this dataframe
  return(pep_data)
}

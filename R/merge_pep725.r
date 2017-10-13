#' Preprocessing of PEP725 data, merges separate files into tidy
#' data, with each observation a line, each column a different
#' parameter value.
#'
#' Mind that the PEP725 database is badly managed. For example the Irish
#' station file includes site description with semi-colons in the description
#' of a semi-colon delimited file format. This will crash the routine. I will
#' not provide a bug fix for these known issues with the underlying
#' PEP725 database due to mismanagement on their part.
#'
#' @param path a path to the PEP725 data (species files only)
#' @return concatted data of all data in the path as a tidy data frame
#' including all normal parameters, and the species name and country code
#' as derived from the file name
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

  print(station_files)

  # merge all the csv files for both the true observations and the
  # station meta-data (name, lat, lon etc.)
  observation_data = do.call(rbind,
                 lapply(data_files, function(file){
                   tmp = utils::read.csv2(file, sep = ";",
                                          stringsAsFactors = FALSE)
                   filename = basename(file)
                   tmp$country = substr(filename,8,9)
                   tmp$species = sub("_"," ",substr(filename,11,nchar(filename)-4))
                   return(tmp)
                   })
                 )

  station_locations = do.call(rbind,
                     lapply(station_files, function(file){
                       print(file)
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
  names(pep_data) = tolower(names(pep_data))

  # return this dataframe
  return(pep_data)
}

#' Preprocessing of PEP725 data, merges separate files into tidy
#' data, with each observation a line, each column a different
#' parameter value.
#'
#' @param path a path to the PEP725 data (either a directory containing tar.gz,
#' unzipped tar.gz files as subdirectories or a single tar.gz file)
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

  # set encoding to circumvent messy encoding in database
  Sys.setlocale("LC_ALL", "C")

  # get tmp directory
  tmpdir = tempdir()

  # if a path to a directory is given, list all tar.gz files
  # if not, assume the linked file is a tar.gz PEP725 file
  # (no formal checks for this are in place)
  if (dir.exists(path.expand(path))){
      archive_files = list.files(path, "^PEP725.*\\.tar\\.gz$",
                                 full.names = TRUE)
  } else {
    # single file evaluation
    if (file.exists(path)){
      archive_files = path
    } else {
      stop("File does not exist, check the filename of the archive file!")
    }
  }

  # return data from a do call, binding the different
  # data sets by row
  do.call("rbind",lapply(archive_files, function(file){

    # unzip all data (thanks Windows for messing up the clean routine)
    untar(file,
          exdir = path.expand(tmpdir))

    # extract only the true data files and station info files
    # drop the BBCH and README data (but don't delete it - delist)
    pep_files = list.files(tmpdir, "PEP725_*", full.names = FALSE)
    pep_files = pep_files[!grepl("^.*PEP725_BBCH.csv$|PEP725_README.txt",
                                 pep_files)]
    data_file = pep_files[!grepl("stations", pep_files)]
    station_file = pep_files[grepl("stations", pep_files)]

    # read in all the required data
    observation_data = utils::read.csv2(sprintf("%s/%s",tmpdir, data_file),
                                        sep = ";",
                                        stringsAsFactors = FALSE)

    station_locations = utils::read.csv2(sprintf("%s/%s",tmpdir,station_file),
                                         sep = ";",
                                         stringsAsFactors = FALSE,
                                         skip = 1,
                                         header = FALSE)

    # discard any columns > 6 (errors in NAME field)
    if(ncol(station_locations) > 6){
      station_locations = station_locations[,-c(7:ncol(station_locations))]
    }

    # manually assign column names to avoid errors with malformed data
    colnames(station_locations) = c("PEP_ID",
                                    "National_ID",
                                    "LON",
                                    "LAT",
                                    "ALT",
                                    "NAME")

    # convert values to numeric and add ancillary data
    station_locations$LON = as.numeric(station_locations$LON)
    station_locations$LAT = as.numeric(station_locations$LAT)

    observation_data$country = substr(data_file, 8, 9)
    observation_data$species = sub("_"," ",substr(data_file,11,nchar(data_file)-4))

    # do a left merge to combine the observational data and the
    # station location meta-data returning basically the original
    # database structure (as a tidy file)
    pep_data = merge(observation_data, station_locations, by = "PEP_ID")
    names(pep_data) = tolower(names(pep_data))

    # cleanup extracted data for good measure
    # and return the combined data frame
    file.remove(list.files(tmpdir, "PEP725_*", full.names = TRUE))
    return(pep_data)
  }))
}

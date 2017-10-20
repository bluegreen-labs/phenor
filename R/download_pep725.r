#' Automatically download PEP725 data
#'
#' Uses the httr and rvest packages to scrape the PEP725 website for
#' phenology data. Login credentials are required so this does not
#' constitute a breach of the license agreement / data policy. Please
#' register here for access to the data:
#' http://www.pep725.eu/data_download/registration.php
#' and carefully read the data use policy:
#' http://www.pep725.eu/downloads/PEP_725_Data_Policy_201012.pdf
#'
#' @param credentials Text file with your PEP725 email and password credentials.
#' The credentials file is a file which lists your email and
#' password on two separate lines (in this order).
#' @param species A species to download, either specified by its
#' species number or species name. list species numbers and names with
#' check_pep725_species(list = TRUE)
#' @param internal completes download internally in a temporary directory and
#' merges the data subsequently using merge_pep725(), returns a nested list
#' of tidy data. internal overrides the path command.
#' @param path the path where to save the downloaded data
#' @return will return csv files of PEP725 data for the selected species
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#' \dontrun{
#' download_pep725(credentials = "~/pep725_login.txt",
#'                 species = 115)
#'}

download_pep725 = function(credentials = NULL,
                           species = 115,
                           path = "~",
                           internal = FALSE){

  # check the validity of the species, return list of
  # numbers to query or stop()
  species_numbers = check_pep725_species(species = species)

  # check if email / password or credential file is available
  if(any(!file.exists(credentials) & is.null(credentials))){
    stop("Credentials file not given or does not exist, check path !")
  } else {
    credentials = as.vector(unlist(read.table(credentials,
                                              stringsAsFactors = FALSE)))
    email = credentials[1]
    password = credentials[2]
  }

  # create login form credentials or generated them from a file
  # the latter is preferred
  login <- list(
    email = email,
    pwd = password,
    submit = "Login"
  )

  # login to set cookie (will not expire until end of session)
  httr::POST("http://www.pep725.eu/login.php", body = login, encode = "form")

  # download data for all species number(s), will merge
  # data if internal = TRUE otherwise return NULL
  all_pep_data = do.call("rbind",lapply(species_numbers, function(number){

    # select the species of interest and pull the table listing
    # all download files
    species_html = httr::POST("http://www.pep725.eu/data_download/data_selection.php",
                              body = list(
                                plant = number,
                                submit1 = "Submit"),
                              encode = "form")

    # extract the links to download
    species_links = xml2::read_html(species_html) %>%
      rvest::html_nodes("td a") %>%
      rvest::html_attr("href")

    # loop over all files for all countries
    # of a particular species, download the (zipped) files
    # select the relevant data, and unzip to specified path
    # use merge_pep725() to merge the downloaded data
    # into a tidy data format
    do.call("rbind",lapply(species_links, function(link){

      # create a temporary file
      tmp = tempfile()

      # download all files for a specific species
      httr::GET(link, httr::write_disk(path = tmp,
                                       overwrite = TRUE))

      # only return data if internal processing is TRUE
      # (file.rename might not be consistent in behaviour
      # hence file.copy() / file.remove() )
      if (internal){
        pep_data = merge_pep725(path = tmp)
        file.remove(tmp)
        return(pep_data)
      } else {
        file.copy(tmp, sprintf("%s/PEP725_%s.tar.gz",
                               path,
                               strsplit(link,"=")[[1]][2]),
                  overwrite = TRUE)
        file.remove(tmp)
      }
    }))
  }))

  # return results if internal flag is set
  # this suppresses file.remove() feedback
  # which
  if (internal){
    return(all_pep_data)
  }
}

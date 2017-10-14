#' Download PEP725 data
#'
#' @param email email used in creating your pep725 download login
#' @param password password as created for the pep725 download login
#' @param credentials file with your email and password credentials, this is
#' the preferred method of storing your credentials over the explicit email and
#' password parameters as these might be forgotten in scripts and publicly
#' uploaded / shared. The credentials file is a file which lists your email and
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
#'
#' \dontrun{
#' download_pep725(credentials = "~/pep725_login.txt,
#'                 species = 115)
#'}

download_pep725 = function(email = NULL,
                           password = NULL,
                           credentials = NULL,
                           species = 115,
                           path = "~",
                           internal = TRUE){

  # check the validity of the species, return list of
  # numbers to query or stop()
  species_numbers = check_pep725_species(species = species)

  # check if email / password or credential file is available
  if(any(is.null(email) | is.null(password)) & is.null(credentials)){
    stop("Matching login credentials (email / password) or a credentials file are required")
  } else {
    credentials = as.vector(unlist(read.table(credentials, stringsAsFactors = FALSE)))
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

  all_pep_data = lapply(species_numbers, function(number){

    # create temporary folder for species processing
    if(internal){
      path = tempdir()
      print(path)
    }

    # select the species of interest and pull the table listing
    # all download files
    species_html = httr::POST("http://www.pep725.eu/data_download/data_selection.php",
                              body = list(
                                plant = number,
                                submit1 = "Submit"),
                              encode = "form")

    # extract the links to download
    species_links = read_html(species_html) %>%
      html_nodes("td a") %>%
      html_attr("href")

    # loop over all files for all countries
    # of a particular species, download the (zipped) files
    # select the relevant data, and unzip to specified path
    # use merge_pep725() to merge the downloaded data
    # into a tidy data format
    lapply(species_links, function(link){

      # create a temporary file
      tmp = tempfile()

      # download all files for a specific species
      httr::GET(link, httr::write_disk(path = tmp, overwrite = TRUE))


    })
  })

  # return all pep data if internal
  # processing is requested
  if(internal){
    return(all_pep_data)
  }
}

test = download_pep725(credentials = "~/pep725login.txt",
                species = c(115,148),
                internal = TRUE)

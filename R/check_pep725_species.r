#' Checks if PEP725 species name or number exists or can be generated
#'
#' @param species A species to download, either specified by its
#' species number or species name.
#' @param list List all species numbers and names as verbose output
#' @return a validated list of species numbers, if not a warning is thrown
#' and any depended routines halted.
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' check_pep725_species(species = 115)
#'}
#'
#' @importFrom magrittr %>%

check_pep725_species <- function(species = NULL,
                                list = FALSE){

  # grab species info from the data selection page
  # this does not require a login and will be used
  # to check the validity of the selected species number
  # or name
  data_selection <- httr::GET("http://www.pep725.eu/data_download/data_selection.php")

  # extract the species numbers and names from the page
  number <- xml2::read_html(data_selection) %>%
    rvest::html_nodes("form select") %>%
    rvest::html_children() %>%
    rvest::html_attr("value") %>%
    as.numeric()

  name <- xml2::read_html(data_selection) %>%
    rvest::html_nodes("form select") %>%
    rvest::html_children() %>%
    rvest::html_text() %>%
    as.character() %>%
    tolower()

  # combine the data in a species info data frame
  species_info <- data.frame(number, name)

  # provide verbose output listing all
  # species names and numbers
  if(list){
    if (is.null(species)){
      return(species_info)
    } else {
      print(species_info,
            row.names = FALSE)
    }
  }

  # if the input is a character vector grep for results
  # this will work on partial matches as well
  if(is.character(species)){
    numbers <- species_info$number[grep(paste(tolower(species),collapse = "|"),
                                    species_info$name)]
    if(length(numbers)==0){
      stop("Species (number) not listed in PEP725 database!")
    } else {
      return(numbers)
    }
  }

  # if the input is a vector of numbers check if they are in
  # the list (and properly numeric)
  if (is.numeric(species) & all(species %in% species_info$number)){
    return(species)
  } else {
    stop("Species (number) not listed in PEP725 database!")
  }
}

check_pep725_species = function(species = 115,
                                list = FALSE){

  # grab species info from the data selection page
  # this does not require a login and will be used
  # to check the validity of the selected species number
  # or name
  data_selection = httr::GET("http://www.pep725.eu/data_download/data_selection.php")

  # extract the species numbers and names from the page
  number = read_html(data_selection) %>%
    html_nodes("form select") %>%
    html_children() %>% html_attr("value") %>% as.numeric()

  name = read_html(data_selection) %>%
    html_nodes("form select") %>%
    html_children() %>% html_text() %>% as.character() %>% tolower()

  # combine the data in a species info data frame
  species_info = data.frame(number,name)

  # provide verbose output listing all
  # species names and numbers
  if(toupper(list)){
    print(species_info, row.names = FALSE)
  }

  # return a complete list of numbers if requested
  # (unlisted feature)
  if ("complete" == tolower(species)){
    return(species_info$number)
  }

  # if the input is a character vector grep for results
  # this will work on partial matches as well
  if(is.character(species)){
    numbers = species_info$number[grep(paste(tolower(species),collapse = "|"),
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

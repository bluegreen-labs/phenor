#' Checks if USA-NPN species exists
#'
#' @param species An USA-NPN species (character or number)
#' @param list List all species numbers and names as verbose output
#' @return a validated list of species numbers, if not a warning is thrown
#' and any depended routines halted.
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' # list all USA-NPN phenophases
#' check_npn_species(species = 3, list = TRUE)
#'}

check_npn_species = function(species = NULL,
                             list = TRUE){
  # download species information
  species_list = jsonlite::fromJSON("http://www.usanpn.org/npn_portal/species/getSpecies.json")

  # drop last column
  species_list = subset(species_list, select = -species)

  # check if the phenophase exists, if provided
  if(!is.null(species)){
    if(any(species %in% species_list$species_id |
           grepl(tolower(species), tolower(paste(species_list$genus, species_list$species)))
    )){

      # print info if requested
      if(list){
        print(species_list[grepl(tolower(species),
                                 tolower(paste(species_list$genus,
                                               species_list$species))) |
                             species_list$species_id %in% species,])
      }
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {

    # if nothing is requested for filtering, just return the full list
    # and a warning
    print(species_list)
    warning("No phenophase provided for further validation.")
  }
}

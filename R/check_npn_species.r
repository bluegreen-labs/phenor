#' Checks if USA-NPN species exists
#'
#' @param species An USA-NPN species (character or number).
#' Will search in both Genus species and common name fields and will match
#' any term within those fields. The search relies on regular expressions so
#' this can be used to be more specific.
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
  species_list <- jsonlite::fromJSON("http://www.usanpn.org/npn_portal/species/getSpecies.json")
  species_list <- subset(species_list, select = -"species_type")

  # check if the genus exists, if provided
  if(!is.null(species)){

    # find species subset results, search (grep) in both Genus species,
    # common name fields and species ID fields
    species_subset <- species_list[species_list$species_id %in% species |
      grepl(tolower(species), tolower(paste(species_list$genus,
                                            species_list$species))) |
      grepl(tolower(species),tolower(species_list$common_name)),]

    # if there is a result return either true or the
    # list of subset values
    if(nrow(species_subset)!= 0){
      if(list){
        return(species_subset)
      } else {
        return(TRUE)
      }
    } else {
      return(FALSE)
    }
  } else {

    # if nothing is requested for filtering, just return the full list
    return(species_list)
  }
}

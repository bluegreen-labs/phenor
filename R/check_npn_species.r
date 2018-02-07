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
  species_list = jsonlite::fromJSON("http://www.usanpn.org/npn_portal/species/getSpecies.json")
  species_list = subset(species_list, select = -species_type)

  # check if the genus exists, if provided
  if(!is.null(species)){

    if(any(species %in% species_list$species_id |
           grepl(tolower(species), tolower(paste(species_list$genus, species_list$species)))
    )){

      # print info if requested
      if(list){
        return(species_list[grepl(tolower(species),
                                 tolower(paste(species_list$genus,
                                               species_list$species))) |
                             species_list$species_id %in% species,])
      } else {
        return(TRUE)
      }
    } else {
      return(FALSE)
    }
  } else {

    # if nothing is requested for filtering, just return the full list
    # and a warning
    return(species_list)
    warning("No phenophase provided for further validation.")
  }
}

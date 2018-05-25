#' Checks if USA-NPN phenophase exists
#'
#' @param phenophase An USA-NPN phenophase (character or number)
#' @param list List all species numbers and names as verbose output
#' @return a validated list of species numbers, if not a warning is thrown
#' and any depended routines halted.
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' # list all USA-NPN phenophases
#' check_npn_phenophases(list = TRUE)
#'}

check_npn_phenophases = function(phenophase = NULL,
                                 list = TRUE){
  # download phenophase information
  phenophases = jsonlite::fromJSON("http://www.usanpn.org/npn_portal/phenophases/getPhenophases.json")

  # check if the phenophase exists, if provided
  if(!is.null(phenophase)){
    if(any(phenophase %in% phenophases$phenophase_id |
           grepl(phenophase, phenophases$phenophase_name))){
      # print info if requested
      if(list){
        print(phenophases[grepl(phenophase, phenophases$phenophase_name) |
                            phenophases$phenophase_id %in% phenophase,])
      }
      return(TRUE)
    } else {
      FALSE
    }
  } else {
    # if nothing is requested for filtering, just return the full list
    # and a warning
    print(phenophases)
    warning("No phenophase provided for further validation.")
  }
}

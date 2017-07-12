#' Null model as defined in
#' Basler et al. 2016 (Agr. For. Meteorlogy)
#' returns the mean across all validation dates
#'
#' @param data: input data (see reference for detailed description),
#' data should be formatted using flat_format()
#' @keywords phenology, model, sequential
#' @export
#' @examples
#'
#' \dontrun{
#' estimate = NLL(data = data)
#'}

null = function(data){
 rep(round(mean(data$transition_dates,na.rm=TRUE)), length(data$transition_dates))
}

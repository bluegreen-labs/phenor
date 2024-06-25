#' Download USA-NPN data using the API interface
#'
#' @param species A species to download, either specified by its
#' species number or species name. list species numbers and names with
#' check_npn_species(list = TRUE)
#' @param phenophase A phenophase number to filter for, a list of phenophases
#' is provided by check_npn_phenophases()
#' @param start start date in YYYY-MM-DD format
#' @param end end date in YYYY-MM-DD format
#' @param extent geographic coordinates constraining the output, defined
#' as bottom left, top right c(lon1, lat1, lon2, lat2) if null returns all
#' data (default = NULL)
#' @param internal completes download internally in a temporary directory
#' @param path the path where to save the downloaded data
#' @return will return a data farme for the selected species,
#' phenophase, temporal and spatial extent or save the data to an RDS file
#' @keywords phenology, model, preprocessing
#' @export
#'
#' @examples
#' \dontrun{
#' pr_dl_npn(species = 3, internal = FALSE)
#'}

pr_dl_npn <- function(
  species = 3,
  phenophase = 371,
  start = "2000-01-01",
  end = "2017-01-01",
  extent = NULL,
  internal = TRUE,
  path = tempdir()
  ){

  # request parameter for logging
  request = "phenor"

  # set url base
  url <- "http://www.usanpn.org/npn_portal/observations/getSummarizedData.json?"

  # formulate the query
  query <- list(species_id = species,
               phenophase_id = phenophase,
               start_date = start,
               end_date = end,
               bottom_left_x1 = extent[1],
               bottom_left_y1 = extent[2],
               upper_rigth_x2 = extent[3],
               upper_right_y2 = extent[4],
               request_src = request)

  # download data using httr
  data <- httr::GET(
    url,
    query = query,
    httr::progress()
    )

  # convert data to a clean data frame
  data <- as.data.frame(jsonlite::fromJSON(httr::content(data, as = "text")))

  # check if the data is valid (contains data)
  # if so return the data frame or write to disk for later processing
  # if not return an error
  if(length(data) == 0){
    stop("Query returned no data, check your input parameters!")
  } else {
    if (!internal){

      # sprintf doesn't deal well with NULL
      if(is.null(phenophase)){
        phenophase <- "NA"
      }

      # write data to file as R data file
      # limiting file size
      saveRDS(data,
              sprintf("%s/phenor_npn_data_%s_%s_%s_%s.rds",
                      path.expand(path),
                      species,
                      phenophase,
                      start,
                      end))
    } else {
      return(data)
    }
  }
}

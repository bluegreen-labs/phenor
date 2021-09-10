#' Download Global CMIP6 driver data
#'
#' Data is downloaded from the Copernicus Data Store, this routine therefore
#' requires a valid account with this service and your ID and key installed on
#' your system using the 'ecmwfr' package.
#'
#' The routine provides access to all available model runs and is a user
#' friendly wrapper around a standard 'ecmwfr' call. Custom data formatting is
#' possible, but herefore we refer to the 'ecmwfr' package.
#'
#' For a full list of the models consult the Copernicus documentation
#' (<https://cds.climate.copernicus.eu/cdsapp#!/dataset/projections-cmip6?tab=overview>) and
#' the data query tools (<https://cds.climate.copernicus.eu/cdsapp#!/dataset/projections-cmip6?tab=form>).
#'
#' @param path a path where to save the gridded data
#' @param end_year end year of the data queried, by default this is set to
#' 2100 but the dataset runs up until 2300
#' @param model CMIP6 model data to download (character vector)
#' @param scenario which RCP scenario to select for (default = "ssp5_8_5",
#' or business as usual). Choose from: historical, ssp1_1_9, ssp1_2.6, ssp4_3_4,
#' ssp5_3_4OS, ssp2_4_5, ssp4_6_0, ssp3_7_0,ssp5_8_5
#' (consult the Copernicus documentation for further details)
#' @param variable which climate variables to download, by default all required
#' variables are downloaded.
#' (default = c("daily_maximum_near_surface_air_temperature",
#' "daily_minimum_near_surface_air_temperature",
#' "precipitation"))
#' @param extent vector with coordinates defining the region of interest defined
#' as ymax, xmin, ymin, xmax in lat/lon (default = c( 90, -180, -90, 180))
#' @param user Copernicus Data Store user ID (a number), on linux do not forget
#' to set options(keyring_backend='file') if you use a file based keyring.
#'
#' @return nothing is returned to the R working environment, files are
#' downloaded and stored on disk
#' @keywords phenology, model, data
#' @export
#' @examples
#'
#' \dontrun{
#' print("bla")
#'}

# create subset of layers to calculate phenology model output on
pr_dl_cmip <- function(
  path = tempdir(),
  end_year = 2100,
  model = "miroc6",
  scenario = "ssp5_8_5",
  variable = c("daily_maximum_near_surface_air_temperature",
               "daily_minimum_near_surface_air_temperature",
               "precipitation"),
  extent = c( 40, -80, 50, -70),
  user
  ){

  if(missing(user)){
    stop("Please provide a valid Copernicus CDS user ID!")
  }

  # check key
  tryCatch(ecmwfr::wf_get_key(
    user = user,
    service = "cds"),
    finally = message(
    "Check your user credentials or your keyring!
    (on linux you might need to set your keyring to a file based system:
    options(keyring_backend='file'))")
    )

  # Downloading data
  lapply(variable, function(var){

    # format request
    request <- list(
      format = "zip",
      temporal_resolution = "daily",
      experiment = scenario,
      level = "single_levels",
      variable = var,
      model = model,
      area = extent,
      date = paste0("2000-01-01/",end_year,"-12-31"),
      dataset_short_name = "projections-cmip6",
      target = "cmip.zip"
    )

    # create phenor temp dir
    dir.create(
      file.path(tempdir(), "phenor"),
      showWarnings = FALSE,
      recursive = TRUE
    )

    ecmwfr::wf_request(
      request = request,
      user = user,
      path = file.path(tempdir(), "phenor")
    )

    # unzipping data in place
    message("unzipping data!")

    unzip(file.path(tempdir(),"phenor", "cmip.zip"),
          exdir = file.path(tempdir(), "phenor"))

    # copy files to final path
    files <- list.files(
      file.path(tempdir(), "phenor"),
      pattern = glob2rx("*.nc"),
      full.names = TRUE)

    file.copy(
      from = files,
      to = path,
      overwrite = TRUE)

    # clean up files in tempdir
    file.remove(
      list.files(file.path(tempdir(), "phenor"),
                 "*",
                 full.names = TRUE)
      )
  })

  # feedback
  message("Download complete!")
}

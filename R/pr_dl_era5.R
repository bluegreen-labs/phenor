#' Download Global ERA5 driver data
#'
#' Data is downloaded from the Copernicus Data Store, this routine therefore
#' requires a valid account with this service and your ID and key installed on
#' your system using the 'ecmwfr' package.
#'
#' The routine provides access to all available model runs and is a user
#' friendly wrapper around a standard 'ecmwfr' call. Custom data formatting is
#' possible, but herefore we refer to the 'ecmwfr' package.
#'
#' @param path a path where to save the gridded data
#' @param year year of the data queried
#' @param extent vector with coordinates defining the region of interest defined
#' as ymax, xmin, ymin, xmax in lat/lon (default = c( 40, -80, 50, -70))
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
pr_dl_era5 <- function(
  path = tempdir(),
  year = 2019,
  extent = c( 40, -80, 50, -70),
  user
  ){

  if(missing(user)){
    stop("Please provide a valid Copernicus CDS user ID!")
  }

  # check key
  tryCatch(ecmwfr::wf_get_key(
    user = as.character(user),
    service = "cds"),
    finally = message(
    "Check your user credentials or your keyring!
    (on linux you might need to set your keyring to a file based system:
    options(keyring_backend='file'))")
    )

  # Downloading data
  request <- list(
    product_type = "reanalysis",
    format = "netcdf",
    variable = c("2m_temperature",
                 "total_precipitation"),
    year = c(year - 1, year),
    month = c("01", "02", "03", "04", "05", "06",
              "07", "08", "09", "10", "11", "12"),
    day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
            "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31"),
    time = c("00:00", "01:00", "02:00", "03:00", "04:00",
             "05:00", "06:00", "07:00", "08:00", "09:00",
             "10:00", "11:00", "12:00", "13:00", "14:00",
             "15:00", "16:00", "17:00", "18:00", "19:00",
             "20:00", "21:00", "22:00", "23:00"),
    area = extent,
    dataset_short_name = "reanalysis-era5-single-levels",
    target = "era5.nc"
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

  # copy files to final path
  files <- list.files(
    file.path(tempdir(), "phenor"),
    pattern = "*\\.nc",
    full.names = TRUE
    )

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

  # feedback
  message("Download complete!")
}

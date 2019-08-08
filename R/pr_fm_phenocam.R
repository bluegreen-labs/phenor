#' Preprocessing of all PhenoCam data into a format which can be ingested
#' by the optimization routines etc. the original nested list is flattened
#' for speed.
#'
#' @param path a path to 1 or 3-day PhenoCam transition date estimates
#' (no validation checks will be done, so mixed files will lead to
#' mixed or duplicated results!)
#' @param direction The phenophase considered, "rising" for spring and
#' "falling" for autumn. When using fall the offset should be 365.
#' (default = rising)
#' @param gcc_value The gcc time series from which the phenophase estimates
#' were derived either (default = gcc_90)
#' @param threshold Threshold (% of seasonal amplitude) used in estimating
#' the phenophases (10, 25, 50) (default = 50)
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param internal return data as structured list to R workspace or write
#' to RDS file (default = TRUE)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' # run with default settings
#' # looks for transition date files derived
#' # through phenocamr in your home directory
#' # change the path to match your setup
#' \dontrun{
#' phenocam_data = pr_fm_phenocam()
#'}

pr_fm_phenocam = function(
  path = tempdir(),
  direction = "rising",
  gcc_value = "gcc_90",
  threshold = 50,
  offset = 264,
  internal = TRUE
  ){

  # helper function to process the data
  format_data <- function(site,
                         transition_files,
                         path,
                         end,
                         metadata){

    # for all sites merge the transition dates if there are multiple files
    # after merging, download the corresponding daymet data and create
    # the parts of the final structured list containing data for further
    # processing
    transition_files_full <- paste(path, transition_files,sep = "/")

    # get individual sites from the filenames
    sites <- unlist(lapply(strsplit(transition_files,"_"),"[[",1))
    files <- transition_files_full[which(sites == site)]

    # merge all transition date data
    data <- do.call("rbind", lapply(files, function(fn){
        data.frame(utils::read.table(fn, header = TRUE, sep = ",") )
      }))

    # if not transition dates exist, return NULL
    if(!any(grepl(threshold,names(data)))){
      return(NULL)
    }

    # throw out all data but the selected gcc value
    data <- data[data$direction == direction &
                  data$gcc_value == gcc_value,
                  grep(threshold,names(data))]

    transition <- as.Date(data[,grep(sprintf("^transition_%s$",threshold),
                                     names(data))])
    lower <- as.Date(data[,grep(sprintf("*%s_lower*",threshold),names(data))])
    upper <- as.Date(data[,grep(sprintf("*%s_upper*",threshold),names(data))])

    # kick out transition dates with large uncertainties (> 30 days)
    # these are most likely false (consider it to be a parameter)
    spread <- abs(lower - upper)
    transition <- transition[spread < 30] # make parameter?

    # grab the location of the site by subsetting the
    site_info <- metadata[which(metadata$site == site),]
    lat <- site_info$lat
    lon <- site_info$lon

    # download daymet data for a given site
    daymet_data <- try(daymetr::download_daymet(
      site = site,
      lat = lat,
      lon = lon,
      start = 1980,
      end = end,
      internal = TRUE,
      silent = TRUE
    )$data,
    silent = TRUE)

    # trap sites outside daymet coverage
    if (inherits(daymet_data,"try-error")){
      return(NULL)
    }

    # calculate the mean daily temperature
    daymet_data$tmean <- (daymet_data$tmax..deg.c. + daymet_data$tmin..deg.c.)/2

    # calculate the long term daily mean temperature
    # and realign it so the first day will be sept 21th (doy 264)
    # and the matching DOY vector
    ltm <- as.vector(by(daymet_data$tmean,
                       INDICES = list(daymet_data$yday),
                       mean,
                       na.rm = TRUE))

    # shift data when offset is < 365
    if (offset < 365){
      ltm = c(ltm[offset:365],ltm[1:(offset - 1)])
      doy_neg = c((offset - 366): -1, 1:(offset - 1))
      doy = c(offset:365,1:(offset - 1))
    } else {
      doy = doy_neg = 1:365
    }

    # how many years to process
    years <- unique(as.numeric(format(transition,"%Y")))

    # create output matrix (holding mean temp.)
    tmean <- matrix(NA,
                   nrow = 365,
                   ncol = length(years))

    # create output matrix (holding min temp.)
    tmin <- matrix(NA,
                  nrow = 365,
                  ncol = length(years))

    # create output matrix (holding max temp.)
    tmax <- matrix(NA,
                  nrow = 365,
                  ncol = length(years))

    # create output matrix (holding vpd)
    vpd <- matrix(NA,
                 nrow = 365,
                 ncol = length(years))

    # create output matrix (holding precip)
    precip <- matrix(NA,
                    nrow = 365,
                    ncol = length(years))

    # create a matrix containing the mean temperature between
    # sept 21th in the previous year until sept 21th in
    # the current year (make this a function parameter)
    # for the default offset, if offset is 365 or larger
    # use the current year only
    for (j in 1:length(years)) {
      if (offset < 365) {
        loc <- which((daymet_data$year == (years[j] - 1) &
                       daymet_data$yday >= offset) |
                     (daymet_data$year == years[j] &
                      daymet_data$yday < offset))
      } else {
        loc <- which(daymet_data$year == years[j])
      }

      tmean[, j] <- daymet_data$tmean[loc]
      tmin[, j] <- daymet_data$tmin..deg.c.[loc]
      tmax[, j] <- daymet_data$tmax..deg.c.[loc]
      precip[, j] <- daymet_data$prcp..mm.day.[loc]
      vpd[, j] <- daymet_data$vp..Pa.[loc]
    }

    # finally select all the transition dates for model validation
    phenophase_years <- as.numeric(format(transition,"%Y"))
    phenophase_doy <- as.numeric(format(transition,"%j"))

    # only select the first instance of a phenophase_doy
    # currently the model frameworks do not handle multiple cycles
    phenophase <- unlist(lapply(years, function(x) {
      phenophase_doy[which(phenophase_years == x)[1]]
    }))

    # calculate daylength
    l <- ncol(tmean)
    Li <- daylength(doy = doy, latitude = lat)
    Li <- matrix(rep(Li,l),length(Li),l)

    # format and return the data
    return(list("site" = site,
                "location" = c(lat,lon),
                "doy" = doy_neg,
                "ltm" = ltm,
                "transition_dates" = phenophase,
                "year" = unique(phenophase_years),
                "Ti" = as.matrix(tmean),
                "Tmini" = as.matrix(tmin),
                "Tmaxi" = as.matrix(tmax),
                "Li" = Li,
                "Pi" = as.matrix(precip),
                "VPDi" = as.matrix(vpd),
                "georeferencing" = NULL
                ))
  }

  # query site list with metadata from the phenocam servers
  metadata <- jsonlite::fromJSON("https://phenocam.sr.unh.edu/webcam/network/siteinfo/")

  # query max year as available through Daymet, lags by a year so
  # subtract 1 year by default. If download fails subtract another year
  end <- as.numeric(format(as.Date(Sys.Date()),"%Y")) - 1

  daymet_test <- try(daymetr::download_daymet(
    start = end,
    end = end,
    internal = TRUE,
    silent = TRUE
  ))

  if (inherits(daymet_test,"try-error")){
    end <- end - 1
  }

  # list all files in the referred path
  transition_files <- list.files(path,
                                "*_transition_dates.csv")

  # get individual sites form the filenames
  sites <- unique(unlist(lapply(strsplit(transition_files,"_"),"[[",1)))

  # track progress
  message(sprintf('Processing %s sites\n', length(sites)))
  pb <- utils::txtProgressBar(min = 0, max = length(sites), style = 3)
  env <- environment()
  i <- 0

  # process data
  validation_data <- lapply(sites, function(x) {
    utils::setTxtProgressBar(pb, i + 1)
    assign("i", i + 1, envir = env)
    format_data(site = x,
                transition_files = transition_files,
                path = path,
                end = end,
                metadata = metadata)
  })

  # close progress bar
  close(pb)

  # rename list variables using the proper site names
  names(validation_data) <- sites

  # assign a class for post-processing
  class(validation_data) <- "phenor_time_series_data"

  # remove out of daymet range sites (prune sites)
  na_loc <- which(is.na(validation_data))
  if (length(na_loc) != 0){
    validation_data = validation_data[-na_loc]
  }

  # return the formatted data
  # either internally or saved as an rds (binary R data file)
  if (internal){
    return(validation_data)
  } else {
    saveRDS(validation_data,
            file = sprintf("%s/phenor_phenocam_data_%s_%s_%s_%s.rds",
                            path,
                            direction,
                            gcc_value,
                            threshold,
                            offset))
  }
}

#' Preprocessing of all PhenoCam data into a format which can be ingested
#' by the optimization routines etc.
#'
#' @param path: a path to 1 or 3-day PhenoCam time series
#' (no validation checks will be done, so mixed files will lead to
#' mixed results!)
#' @param direction: rising = spring, falling = autumn
#' @param gcc_valuel: gcc_90, gcc_mean, gcc_50 etc.
#' @param transition: 10, 25, 50 = default (threshold)
#' @param offset: offset of the time series in DOY (default = 264, sept 21)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples

setwd("~/Dropbox/Data/tmp/transition_dates/")

process.phenocam = function(path = ".",
                            direction = "rising",
                            gcc_value = "gcc_90",
                            transition = 50,
                            offset = 264){

  # helper function to process the data
  format_data = function(site, transition_files, metadata){

    # for all sites merge the transition dates if there are multiple files
    # after merging, download the corresponding daymet data and create
    # the parts of the final structured list containing data for further
    # processing

    # get individual sites form the filenames
    sites = unlist(lapply(strsplit(transition_files,"_"),"[[",1))
    files = transition_files[which(sites == site)]

    # merge all transition date data
    data = do.call("rbind", lapply(files, function(fn)
      data.frame( read.table(fn, header = TRUE, sep = ",") )))

    # throw out all data but gcc_90
    data = data[data$direction == direction &
                  data$gcc_value == gcc_value,
                  grep(transition,names(data))]

    # kick out transition dates with large uncertainties (> 30 days)
    # these are most likely false (consider it to be a parameter)
    spread = abs(as.Date(data$transition_50_lower_ci) - as.Date(data$transition_50_upper_ci))
    data = data[spread < 30,]

    # grab the location of the site by subsetting the
    site_info = metadata[which(metadata$site == site),]
    lat = site_info$lat
    lon = site_info$lon

    # min and max range of the phenology data
    # -1 for min_year as we need data from the previous year for cold
    # hardening
    start_yr = as.numeric(min(format(as.Date(data$transition_50),"%Y"))) - 1
    end_yr = as.numeric(max(format(as.Date(data$transition_50),"%Y")))

    # download daymet data for a given site
    daymet_data = download.daymet(
      site = i,
      lat = lat,
      lon = lon,
      start_yr = 1980,
      end_yr = end_yr,
      internal = "data.frame",
      quiet = TRUE
    )$data

    # calculate the mean daily temperature
    daymet_data$tmean = (daymet_data$tmax..deg.c. + daymet_data$tmin..deg.c.)/2

    # calculate the long term daily mean temperature and realign it so the first
    # day will be sept 21th (doy 264) and the matching DOY vector
    ltm = as.vector(by(daymet_data$tmean, INDICES = list(daymet_data$yday), mean))
    ltm = c(ltm[offset:365],ltm[1:(offset - 1)])
    doy = c(offset:365,1:(offset - 1))

    # slice and dice the data
    years = as.numeric(format(as.Date(data$transition_50),"%Y"))

    # create output matrix (holding temperature)
    temperature = matrix(NA,
                         nrow = 365,
                         ncol = length(years))

    # create a matrix containing the mean temperature between
    # sept 21th in the previous year until sept 21th in
    # the current year (make this a function parameter)
    for (j in 1:length(years)) {
      temperature[,j] = subset(daymet_data, (year == (years[j] - 1) & yday >= offset)|
                          ( year == years[j] & yday < offset ) )$tmean
    }

    # finally select all the transition dates for model validation
    phenophase = as.numeric(format(as.Date(data$transition_50),"%j"))

    # format the data
    data = list("location" = c(lat,lon),
                "doy" = doy,
                "ltm" = ltm,
                "phenophase" = phenophase,
                "tmean" = temperature)

    # return the formatted data
    return(data)
  }

  # query site list with metadata from the phenocam servers
  metadata = jsonlite::fromJSON("https://phenocam.sr.unh.edu/webcam/network/siteinfo/")

  # list all files in the referred path
  transition_files = list.files(path,"*_transition_dates.csv")

  # get individual sites form the filenames
  sites = unique(unlist(lapply(strsplit(transition_files,"_"),"[[",1)))

  # construct validation data using the helper function
  # format_data() above
  validation_data = lapply(sites, function(x) {
    format_data(site = x,
                transition_files = transition_files,
                metadata = metadata)
  })

  # rename list variables using the proper site names
  names(validation_data) = sites

  # return the formatted data
  return(validation_data)
}

bla = process.phenocam()
print(str(bla))
saveRDS(bla, "/data/Dropbox/phenor_validation_data.rds")

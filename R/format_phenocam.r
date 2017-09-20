#' Preprocessing of all PhenoCam data into a format which can be ingested
#' by the optimization routines etc. the original nested list is flattened
#' for speed.
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
#'
#' \dontrun{
#' # run with default settings
#' # looks for transition date files derived
#' # through phenocamr in your home directory
#' # change the path to match your setup
#' phenocam_data = format_phenocam()
#'}

format_phenocam = function(path = "~",
                           direction = "rising",
                           gcc_value = "gcc_90",
                           threshold = 50,
                           offset = 264){

  # helper function to process the data
  format_data = function(site,
                         transition_files,
                         path,
                         end,
                         metadata){

    # provide some feedback
    cat(sprintf('Processing site: %s !\n', site))

    # for all sites merge the transition dates if there are multiple files
    # after merging, download the corresponding daymet data and create
    # the parts of the final structured list containing data for further
    # processing
    transition_files_full = paste(path, transition_files,sep = "/")

    # get individual sites form the filenames
    sites = unlist(lapply(strsplit(transition_files,"_"),"[[",1))
    files = transition_files_full[which(sites == site)]

    # merge all transition date data
    data = do.call("rbind", lapply(files, function(fn)
      data.frame( read.table(fn, header = TRUE, sep = ",") )))

    if(!any(grepl(threshold,names(data)))){
      cat(sprintf('No transition dates for threshold %s
                  at site %s !\n',threshold, site))
      return(NA)
    }

    # throw out all data but the selected gcc value
    data = data[data$direction == direction &
                  data$gcc_value == gcc_value,
                  grep(threshold,names(data))]

    transition = as.Date(data[,grep(sprintf("^transition_%s$",threshold),names(data))])
    lower = as.Date(data[,grep(sprintf("*%s_lower*",threshold),names(data))])
    upper = as.Date(data[,grep(sprintf("*%s_upper*",threshold),names(data))])

    # kick out transition dates with large uncertainties (> 30 days)
    # these are most likely false (consider it to be a parameter)
    spread = abs(lower - upper)
    transition = transition[spread < 30] # make parameter?

    # grab the location of the site by subsetting the
    site_info = metadata[which(metadata$site == site),]
    lat = site_info$lat
    lon = site_info$lon

    # download daymet data for a given site
    daymet_data = try(daymetr::download_daymet(
      site = site,
      lat = lat,
      lon = lon,
      start = 1980,
      end = end,
      internal = "data.frame",
      quiet = TRUE
    )$data)

    # trap sites outside daymet coverage
    if (inherits(daymet_data,"try-error")){
      cat(sprintf('Site: %s is located outside Daymet coverage
                      will be pruned!\n', site))
      return(NA)
    }

    # calculate the mean daily temperature
    daymet_data$tmean = (daymet_data$tmax..deg.c. + daymet_data$tmin..deg.c.)/2

    # calculate the long term daily mean temperature
    # and realign it so the first day will be sept 21th (doy 264)
    # and the matching DOY vector
    ltm = as.vector(by(daymet_data$tmean,
                       INDICES = list(daymet_data$yday),
                       mean,
                       na.rm = TRUE))

    # shift data when offset is < 365
    if (offset < 365){
      ltm = c(ltm[offset:365],ltm[1:(offset - 1)])
      doy_neg = c((offset - 366):-1,1:(offset - 1))
      doy = c(offset:365,1:(offset - 1))
    } else {
      doy = doy_neg = 1:365
    }

    # slice and dice the data
    years = unique(as.numeric(format(transition,"%Y")))

    # create output matrix (holding mean temp.)
    tmean = matrix(NA,
                   nrow = 365,
                   ncol = length(years))

    # create output matrix (holding min temp.)
    tmin = matrix(NA,
                  nrow = 365,
                  ncol = length(years))

    # create output matrix (holding max temp.)
    tmax = matrix(NA,
                  nrow = 365,
                  ncol = length(years))

    # create output matrix (holding vpd)
    vpd = matrix(NA,
                 nrow = 365,
                 ncol = length(years))

    # create output matrix (holding precip)
    precip = matrix(NA,
                    nrow = 365,
                    ncol = length(years))

    # create a matrix containing the mean temperature between
    # sept 21th in the previous year until sept 21th in
    # the current year (make this a function parameter)
    # for the default offset, if offset is 365 or larger
    # use the current year only
    for (j in 1:length(years)) {
      if (offset < 365) {
        tmean[, j] = subset(daymet_data,
                                  (year == (years[j] - 1) &
                                     yday >= offset) |
                                    (year == years[j] &
                                       yday < offset))$tmean

        tmin[, j] = subset(daymet_data,
                                  (year == (years[j] - 1) &
                                     yday >= offset) |
                                    (year == years[j] &
                                       yday < offset))$tmin..deg.c.

        tmax[, j] = subset(daymet_data,
                           (year == (years[j] - 1) &
                              yday >= offset) |
                             (year == years[j] &
                                yday < offset))$tmax..deg.c.

        precip[, j] = subset(daymet_data,
                           (year == (years[j] - 1) & yday >= offset) |
                             (year == years[j] &
                                yday < offset))$prcp..mm.day.
        vpd[, j] = subset(daymet_data,
                             (year == (years[j] - 1) & yday >= offset) |
                               (year == years[j] &
                                  yday < offset))$vp..Pa.

      } else {
        tmean[, j] = subset(daymet_data, year == years[j])$tmean
        tmin[, j] = subset(daymet_data, year == years[j])$tmin..deg.c.
        tmax[, j] = subset(daymet_data, year == years[j])$tmax..deg.c.
        precip[, j] = subset(daymet_data, year == years[j])$prcp..mm.day.
        vpd[, j] = subset(daymet_data, year == years[j])$vp..Pa.
      }
    }

    # finally select all the transition dates for model validation
    phenophase_years = as.numeric(format(transition,"%Y"))
    phenophase_doy = as.numeric(format(transition,"%j"))

    # only select the first instance of a phenophase_doy
    # currently the model frameworks do not handle multiple cycles
    phenophase = unlist(lapply(years, function(x) {
      phenophase_doy[which(phenophase_years == x)[1]]
    }))

    # calculate daylength
    l = ncol(tmean)
    Li = daylength(doy = doy, latitude = lat)
    Li = matrix(rep(Li,l),length(Li),l)

    # format the data
    data = list("site" = site,
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
                "VPDi" = as.matrix(vpd)
                )

    # assign a class for post-processing
    class(data) = "phenor_time_series_data"

    # return the formatted data
    return(data)
  }

  # query site list with metadata from the phenocam servers
  metadata = jsonlite::fromJSON("https://phenocam.sr.unh.edu/webcam/network/siteinfo/")

  # query max year as available through Daymet, lags by a year so
  # subtract 1 year by default. If download fails subtract another year
  end = as.numeric(format(as.Date(Sys.Date()),"%Y")) - 1

  daymet_test = try(daymetr::download_daymet(
    start = end,
    end = end,
    internal = "data.frame",
    quiet = TRUE
  ))

  if (inherits(daymet_test,"try-error")){
    end = end - 1
  }

  # list all files in the referred path
  transition_files = list.files(path,
                                "*_transition_dates.csv")

  # get individual sites form the filenames
  sites = unique(unlist(lapply(strsplit(transition_files,"_"),"[[",1)))

  # construct validation data using the helper function
  # format_data() above
  validation_data = lapply(sites, function(x) {
    format_data(site = x,
                transition_files = transition_files,
                path = path,
                end = end,
                metadata = metadata)
  })

  # rename list variables using the proper site names
  names(validation_data) = sites

  # remove out of daymet range sites (prune sites)
  na_loc = which(is.na(validation_data))
  if (length(na_loc) != 0){
    validation_data = validation_data[-na_loc]
  }

  # return the formatted data
  return(validation_data)
}

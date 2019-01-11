#' Pre-processing of USA-NPN data
#'
#' Combines data into a format which can be ingested
#' by the optimization routines. Data is aggregated on an
#' individual_id basis, rather than on the user defined site_id.
#'
#' @param data an USA-NPN data frame as returned by download_npn()
#' @param phenophase a phenophase to include as validation statistic
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param internal return data as structured list to R workspace or write
#' to RDS file (default = TRUE)
#' @param path path where to save the generated data file
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' npn_data = format_npn()
#'}

format_npn = function(
  data,
  phenophase = 371,
  offset = 264,
  path = tempdir(),
  internal = TRUE
  ){

  # check
  if (missing(data)){
    stop("No data provided")
  }

  # helper function to process the data
  format_data = function(site){

    # subset the original data
    data_subset = data[data$individual_id == site &
                       data$phenophase_id == phenophase,]

    # if nothing in subset return NULL
    if(nrow(data_subset)==0){
      return(NULL)
    }

    # grab the location of the site by subsetting the
    lat = unique(data_subset$latitude)
    lon = unique(data_subset$longitude)


    # get mean transition dates per year by site
    mean_doy = by(data_subset$first_yes_doy,
                  INDICES = data_subset$first_yes_year,
                  mean,
                  na.rm = TRUE)
    transition_dates = round(as.vector(mean_doy))

    # get years matching transition dates
    years = as.numeric(names(mean_doy))

    # download daymet data for a given site
    daymet_data = try(daymetr::download_daymet(
      site = site,
      lat = lat,
      lon = lon,
      start = 1980,
      end = end,
      internal = TRUE,
      silent = TRUE
    )$data)

    # trap sites outside daymet coverage
    if (inherits(daymet_data,"try-error")){
      #cat(sprintf('Site: %s outside Daymet coverage will be pruned!\r\r', site))
      return(NULL)
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
                            ("year" == (years[j] - 1) & "yday" >= offset) |
                              ("year" == years[j] &
                                 "yday" < offset))$tmean

        tmin[, j] = subset(daymet_data,
                           ("year" == (years[j] - 1) & "yday" >= offset) |
                             ("year" == years[j] &
                                "yday" < offset))$tmin..deg.c.

        tmax[, j] = subset(daymet_data,
                           ("year" == (years[j] - 1) & "yday" >= offset) |
                             ("year" == years[j] &
                                "yday" < offset))$tmax..deg.c.

        precip[, j] = subset(daymet_data,
                             ("year" == (years[j] - 1) & "yday" >= offset) |
                               ("year" == years[j] &
                                  "yday" < offset))$prcp..mm.day.
        vpd[, j] = subset(daymet_data,
                          ("year" == (years[j] - 1) & "yday" >= offset) |
                            ("year" == years[j] &
                               "yday" < offset))$vp..Pa.
      } else {
        tmean[, j] = subset(daymet_data, "year" == years[j])$tmean
        tmin[, j] = subset(daymet_data, "year" == years[j])$tmin..deg.c.
        tmax[, j] = subset(daymet_data, "year" == years[j])$tmax..deg.c.
        precip[, j] = subset(daymet_data, "year" == years[j])$prcp..mm.day.
        vpd[, j] = subset(daymet_data, "year" == years[j])$vp..Pa.
      }
    }

    # calculate daylength
    l = ncol(tmean)
    Li = daylength(doy = doy, latitude = lat)
    Li = matrix(rep(Li,l),length(Li),l)

    # format and return the data
    return(list("site" = as.character(site),
                "location" = c(lat,lon),
                "doy" = doy_neg,
                "ltm" = ltm,
                "transition_dates" = transition_dates,
                "year" = years,
                "Ti" = as.matrix(tmean),
                "Tmini" = as.matrix(tmin),
                "Tmaxi" = as.matrix(tmax),
                "Li" = Li,
                "Pi" = as.matrix(precip),
                "VPDi" = as.matrix(vpd),
                "georeferencing" = NULL
    ))
  }

  # query max year as available through Daymet, lags by a year so
  # subtract 1 year by default. If download fails subtract another year
  end = as.numeric(format(as.Date(Sys.Date()),"%Y")) - 1

  daymet_test = try(daymetr::download_daymet(
    start = end,
    end = end,
    internal = TRUE,
    silent = TRUE
  ))

  if (inherits(daymet_test,"try-error")){
    end = end - 1
  }

  # get individual sites form the filenames
  # a "site" in this case is defined as a single
  # individual (or location) for which recurrent data is observed
  sites = unique(data$individual_id)

  # track progress
  pb = utils::txtProgressBar(min = 0, max = length(sites), style = 3)
  env = environment()
  i = 0
  cat(sprintf('Processing %s individuals (sites)\n', length(sites)))

  # process data
  validation_data = lapply(sites, function(x) {
    utils::setTxtProgressBar(pb, i + 1)
    assign("i", i+1, envir = env)
    format_data(site = x)
  })

  # close progress bar
  close(pb)

  # rename list variables using the proper site names
  names(validation_data) = sites

  # assign a class for post-processing
  class(validation_data) = "phenor_time_series_data"

  # remove out of daymet range sites (prune sites)
  na_loc = which(is.na(validation_data))
  if (length(na_loc) != 0){
    validation_data = validation_data[-na_loc]
  }

  # return the formatted data
  # either internally or saved as an rds (binary R data file)
  if (internal){
    return(validation_data)
  } else {
    saveRDS(validation_data,
            file = sprintf("%s/phenor_npn_data_%s_%s.rds",
                           path,
                           phenophase,
                           offset))
  }
}

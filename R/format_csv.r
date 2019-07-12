#' Pre-processing of CSV based data
#'
#' Combines CSV data into a format which can be ingested
#' by the optimization routines. The CSV requires columns named:
#' ID, lat, lon, phenophase, year, doy
#'
#' subsets of the phenophases can be made, or are ignored if
#' not specified.
#'
#' @param file an CSV file with phenology observation dates
#' @param phenophase a phenophase to include as validation statistic
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param internal return data as structured list to R workspace or write
#' to RDS file (default = TRUE)
#' @param path output directory for converted data
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' npn_data = format_csv()
#'}

format_csv <- function(
  file = "~/Desktop/your_phenology.csv",
  phenophase,
  offset = 264,
  internal = TRUE,
  path = tempdir()
  ){

  # helper function to process the data
  format_data = function(site){

    # throw out all data but the selected site
    df = data[data$ID == site,]

    #check for unique
    if (length(unique(df$lat))>1){
      stop ("individual sites span more than one latitude")
    }
    if (length(unique(df$long))>1){
      stop ("individual sites span more than one longitude")
    }

    # grab the location of the site
    lat = unique(df$lat)
    lon = unique(df$long)

    # download daymet data for a given site
    daymet_data = try(daymetr::download_daymet(
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
    years = unique(df$year)

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

      # I need the preceding year for modelling
      # purposes the first usable year is therefore 1981
      if(years[j] == 1980){
        return(NULL)
      }

      if (offset < 365) {
        tmean[, j] = subset(daymet_data,
                            (year == (years[j] - 1) & yday >= offset) |
                              (year == years[j] &
                                 yday < offset))$tmean

        tmin[, j] = subset(daymet_data,
                            (year == (years[j] - 1) & yday >= offset) |
                              (year == years[j] &
                                 yday < offset))$tmin..deg.c.

        tmax[, j] = subset(daymet_data,
                           (year == (years[j] - 1) & yday >= offset) |
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
    phenophase_years = df$year
    phenophase_doy = df$doy

    # only select the first instance of a phenophase_doy
    # currently the model frameworks do not handle multiple cycles
    phenophase_obs = unlist(lapply(years, function(x) {
      phenophase_doy[which(phenophase_years == x)[1]]
    }))

    # calculate daylength
    l = ncol(tmean)
    Li = daylength(doy = doy, latitude = lat)
    Li = matrix(rep(Li,l),length(Li),l)

    # format and return the data
    return(list("site" = site,
                "location" = c(lat,lon),
                "doy" = doy_neg,
                "ltm" = ltm,
                "transition_dates" = phenophase_obs,
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

  # read in the data from csv and subset by phenological stage
  data = utils::read.table(file,
                    sep = ",",
                    header = TRUE,
                    stringsAsFactors = FALSE)

  #fail if phenophase not contained in data
  if(!phenophase%in%data$phenophase){
    stop(paste0("dataset does not contain any observations of the phenophase ", phenophase))
  }

  # subset if a phenophase is specified
  if(missing(phenophase)){
    stop("please specify a phenophase to process")
  }else{
    data=data[data$phenophase==phenophase,]
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
  sites = unique(data$ID)

  # track progress
  cat(sprintf('Processing %s sites\n', length(sites)))
  pb = utils::txtProgressBar(min = 0, max = length(sites), style = 3)
  env = environment()
  i = 0

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
            file = sprintf("%s/phenor_vame_data_%s_%s.rds",
                           path,
                           phenophase,
                           offset))
  }
}

#' Preprocessing of all PhenoCam data into a format which can be ingested
#' by the optimization routines etc. Ingests MODIS phenology dates as ingested
#' through the python based Google Earth Engine (GEE) subset tool. This is a
#' fallback method in the absence of a working version of R MODISTools.
#' (https://khufkens.github.io/gee_subset/)
#'
#' @param path a path to GEE provided MCD12Q2 phenology dates
#' @param phenophase Phenological phase, Increase(1), Maximum,
#' Decrease or Minimum (default = Increase)
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' modis_data = format_gee()
#'}

format_gee = function(path = "~",
                         phenophase = "Increase",
                         offset = 264){

  # helper function to process the data
  format_data = function(site, transition_files, path){

    # define the start date of the MODIS
    # phenology product (days since Jan. 1th 2000)
    start_date = as.Date("2000-01-01")

    # get individual sites form the filenames
    sites = unlist(lapply(strsplit(basename(transition_files),"_"),"[[",1))
    file = transition_files[which(sites == site)]

    # merge all transition date data
    modis_data = utils::read.table(file, header = TRUE, sep = ",")

    # grab the site years from the product name
    years = as.numeric(format(as.Date(modis_data$date),"%Y"))

    # grab the location of the site
    lat = unique(modis_data$latitude)
    lon = unique(modis_data$longitude)

    # grab the values, set NA flag and take the median
    modis_data = modis_data[,grep(phenophase,colnames(modis_data))]
    modis_data = unlist(by(modis_data, INDICES = list(years), stats::median, na.rm = TRUE))

    # count the days from the start date, and report the DOY
    phenophase = as.numeric(format(start_date + modis_data, "%j"))

    # min and max range of the daymet data to pull
    # -1 for min_year as we need data from the previous year for cold
    # hardening
    start = min(years) - 1
    end = max(years)

    # download daymet data for a given site
    daymet_data = try(daymetr::download_daymet(
      site = site,
      lat = lat,
      lon = lon,
      start = 1980,
      end = end,
      internal = TRUE,
      quiet = TRUE
    )$data,
    silent = TRUE)

    # trap sites outside daymet coverage
    if (inherits(daymet_data,"try-error")){
      return(NULL)
    }

    # calculate the mean daily temperature
    daymet_data$tmean = (daymet_data$tmax..deg.c. + daymet_data$tmin..deg.c.)/2

    # calculate the long term daily mean temperature and realign it so the first
    # day will be sept 21th (doy 264) and the matching DOY vector
    ltm = as.vector(by(daymet_data$tmean, INDICES = list(daymet_data$yday), mean))
    ltm = c(ltm[offset:365],ltm[1:(offset - 1)])

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

    # calculate daylength
    l = ncol(tmean)
    Li = daylength(doy = doy, latitude = lat)
    Li = matrix(rep(Li,l),length(Li),l)

    # recreate the validation data structure (new format)
    # but with concatted data
    data = list("site" = site,
                "location" = c(lat,lon),
                "doy" = doy_neg,
                "ltm" = ltm,
                "transition_dates" = phenophase,
                "year" = years,
                "Ti" = as.matrix(tmean),
                "Tmini" = as.matrix(tmin),
                "Tmaxi" = as.matrix(tmax),
                "Li" = Li,
                "Pi" = as.matrix(precip),
                "VPDi" = as.matrix(vpd),
                "georeferencing" = NULL
    )

    # assign a class for post-processing
    class(data) = "phenor_time_series_data"

    # return the formatted data
    return(data)
  } # end of helper function

  #----

  # list all files in the referred path
  transition_files = list.files(path,"*_MCD12Q2_gee_subset.csv",
                                full.names = TRUE)

  # get individual sites form the filenames
  sites = unique(unlist(lapply(strsplit(basename(transition_files),"_"),"[[",1)))

  # track progress
  cat(sprintf('Processing %s sites\n', length(sites)))
  pb = txtProgressBar(min = 0, max = length(sites), style = 3)
  env = environment()
  i = 0

  # get data
  validation_data = lapply(sites, function(x) {
    setTxtProgressBar(pb, i + 1)
    assign("i", i+1, envir = env)
    format_data(site = x,
                transition_files = transition_files,
                path = normalizePath(path))
  })

  # rename list variables using the proper site names
  names(validation_data) = sites

  # assign a class for post-processing
  class(validation_data) = "phenor_time_series_data"

  # return the formatted data
  return(validation_data)
}

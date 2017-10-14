#' Preprocessing of all PhenoCam data into a format which can be ingested
#' by the optimization routines etc.
#'
#' @param path a path to MODISTools MCD12Q2 phenology dates
#' @param direction Increase, Maximum, Decrease or Minimum (default = Increase)
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' modis_data = format_modis()
#'}

process_modis = function(path = "~",
                         direction = "Increase",
                         offset = 264){

  # helper function to process the data
  format_data = function(site, transition_files, path){

    # for all sites merge the transition dates if there are multiple files
    # after merging, download the corresponding daymet data and create
    # the parts of the final structured list containing data for further
    # processing

    # get individual sites form the filenames
    sites = unlist(lapply(strsplit(transition_files,"_"),"[[",1))
    transition_files_full = paste(path, transition_files,sep = "/")
    files = transition_files_full[which(sites == site)]

    # merge all transition date data
    modis_data = utils::read.table(files, header = FALSE, sep = ",")

    # grab the site years from the product name
    years = unique(as.numeric(substring(modis_data[,8],2,5)))

    # grab the location of the site
    lat_lon = unlist(strsplit(as.character(modis_data[1,9]),"Lon"))
    lat = as.numeric(gsub("Lat","",lat_lon[1]))
    lon = as.numeric(lapply(strsplit(lat_lon[2],"Samp"),"[[",1))

    # throw out all data but gcc_90
    modis_data = modis_data[grep(direction, as.character(modis_data[,6])),
                            11:ncol(modis_data)]
    modis_data[modis_data == 32767] = NA
    modis_data = round(apply(modis_data, 1, stats::median, na.rm = TRUE))

    # min and max range of the phenology data
    # -1 for min_year as we need data from the previous year for cold
    # hardening
    start = min(years) - 1
    end = max(years)

    # download daymet data for a given site
    daymet_data = daymetr::download_daymet(
      site = site,
      lat = lat,
      lon = lon,
      start = 1980,
      end = end_yr,
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
    modis_doy = as.numeric(format(seq(as.Date("2001/1/1"),Sys.Date(), "days"),"%j"))
    phenophase = modis_doy[modis_data]

  # recreate the validation data structure (new format)
  # but with concatted data
  data = list("site" = site,
              "location" = c(lat, lon),
              "doy" = doy_neg,
              "ltm" = ltm,
              "transition_dates" = phenophase,
              "year" = pep_subset$year,
              "Ti" = Ti, # temperature
              "Tmini" = Tmini,
              "Tmaxi" = Tmaxi,
              "Li" = Li,
              "Pi" = Pi,
              "VPDi" = VPDi,
              "georeferencing" = NULL
  )

    # assign a class for post-processing
    class(data) = "phenor_time_series_data"

    # return the formatted data
    return(data)
  }

  # list all files in the referred path
  transition_files = list.files(path,"*_MCD12Q2.asc")

  # get individual sites form the filenames
  sites = unique(unlist(lapply(strsplit(transition_files,"_"),"[[",1)))

  # construct validation data using the helper function
  # format_data() above
  validation_data = lapply(sites, function(x) {
    format_data(site = x,
                transition_files = transition_files,
                path = path)
  })

  # rename list variables using the proper site names
  names(validation_data) = sites

  # assign a class for post-processing
  class(validation_data) = "phenor_time_series_data"

  # return the formatted data
  return(validation_data)
}

#' Preprocessing of PEP725 data into a format which can be ingested
#' by the optimization routines, the original nested list is flattened
#' for speed.
#'
#' Some pre-processing steps are required as downloading the PEP725 data is a
#' mess (this database needs an API). So, for the species of interest
#' download the separate zipped files. Unzip all files and put the respective
#' scientific data in one folder.
#'
#' @param path: a path to the PEP725 data (species files only)
#' @param bbch: which bbch
#' @param species: species to select from merged PEP725 file
#' @param offset: offset of the time series in DOY (default = 264, sept 21)
#' @param count: minimum number of acquisitions per location
#' @return returns a nested list of site locations, their respective
#' phenological metrics and matching environmental data as extracted from
#' the E-OBS product (corrected for altitude using a lapse rate of 5C/km.)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' \dontrun{
#' # run with default settings
#' # looks for transition date files derived
#' # through phenocamr in your home directory
#' # change the path to match your setup
#' phenocam_data = format_pep725()
#'}

format_pep725 = function(pep_path = "~",
                         eobs_path = "~",
                         bbch = "11",
                         species = NULL,
                         offset = 264,
                         count = 60,
                         resolution = 0.25){

  # helper function to format the data
  # for a given site
  format_data = function(site = site,
                         offset = offset){

    # subset the data based upon the year to evaluate
    # necessary for E-OBS data subsetting
    pep_subset = pep_data[which(pep_data$pep_id == site),]

    # create spatial object with
    points = sp::SpatialPoints(cbind(pep_subset$lon[1],
                                     pep_subset$lat[1]),
                  proj4string = sp::CRS(raster::projection(eobs_data[[1]])))

    # subset raster data (only once, then replicate data for speed)
    # [subsetting the raster data is a bottleneck]
    temperature = raster::extract(eobs_data[[1]], points)

    # close to water bodies raster cells can be empty and return
    # NA, in this case set the output to NULL and remove these
    # values later on
    if (all(is.na(temperature))){
      return(NULL)
    } else {
      precipitation = raster::extract(eobs_data[[2]], points)
      temperature = matrix(rep(temperature, nrow(pep_subset)),
                           length(temperature),
                           nrow(pep_subset))
      precipitation = matrix(rep(precipitation, nrow(pep_subset)),
                             length(precipitation),
                             nrow(pep_subset))
    }

    # calculate long term mean temperature
    ltm = as.vector(unlist(by(temperature[which(years >= 1980),1],
             INDICES = yday[which(years >= 1980)],
             mean,
             na.rm = TRUE)))[1:365]

    # calculate lapse rate and correct temperatures
    lapse_rate = as.vector((raster::extract(eobs_data[[3]], points[1]) - pep_subset$alt[1]) * 0.005)
    temperature = rbind(temperature + lapse_rate, pep_subset$year)
    ltm = ltm + lapse_rate

    # now for all years create subsets for temperature
    # and precipitation
    Ti = apply(temperature, 2, function(x){
      layers = which((years == (x[length(x)] - 1) & yday >= offset) |
                       (years == x[length(x)] & yday < offset))[1:365]
      return(x[layers])
    })

    Pi = apply(rbind(precipitation, pep_subset$year), 2, function(x){
      layers = which((years == (x[length(x)] - 1) & yday >= offset) |
                       (years == x[length(x)] & yday < offset))[1:365]
      return(x[layers])
    })

    # create a doy vector
    doy = c(offset:365,1:(offset - 1))

    # calculate daylength
    l = ncol(Ti)
    Li = daylength(doy = doy, latitude = pep_subset$lat[1])
    Li = matrix(rep(Li,l),length(Li),l)

    # recreate the validation data structure (new format)
    # but with concatted data
    data = list("site" = site,
                "location" = c(pep_subset$lat[1], pep_subset$lon[1]),
                "doy" = doy,
                "ltm" = ltm,
                "transition_dates" = pep_subset$day,
                "year" = pep_subset$year,
                "Ti" = Ti,
                "Tmini" = NULL,
                "Tmaxi" = NULL,
                "Li" = Li,
                "Pi" = Pi,
                "VPDi" = NULL
    )

    # assign a class for post-processing
    class(data) = "phenor_time_series_data"

    # return the list
    return(data)
  }

  cat("* Merging PEP725 data files in: \n")
  cat(sprintf("  %s\n", pep_path))
  # concat the raw PEP725 data
  pep_data = merge_pep725(path = pep_path)

  # removing out of E-OBS climate data range
  # PEP725 observations
  cat(" |_ removing data out of range of the E-OBS climate data \n")
  pep_data = pep_data[which(pep_data$year >= 1950 &
                            pep_data$year <= max(pep_data$year)),]

  # filtering on species
  if (is.null(species)){
    cat(" |_ including all species \n")
  } else {
    cat(sprintf(" |_ selecting species: %s \n", species))
    pep_data = pep_data[which(pep_data$species == species),]
  }

  # filtering based upon bbch value (phenophase)
  cat(sprintf(" |_ selecting phenophase: %s \n", bbch))
  pep_data = pep_data[which(pep_data$bbch == bbch),]

  cat(sprintf(" |_ excluding sites with less than: %s observations \n", count))
  # count number of observations for each
  # unique pep725 site
  sites = unique(pep_data$pep_id)
  counts = unlist(lapply(sites,function(x){
      length(which(pep_data$pep_id == x))
      })
    )

  # make a selection based upon minimum number
  # of observations required per site
  selection = sites[which(counts >= count)]
  pep_data = pep_data[pep_data$pep_id %in% selection,]

  # sanity check on what remains
  if (nrow(pep_data) == 0 ){
    stop("no data remaining after screening for the requested criteria
         check your species name and observation count restrictions!")
  }

  # get remaining unique years and sites
  sites = unique(pep_data$pep_id)
  years = unique(pep_data$year)

  # loading E-OBS data for subsetting
  cat(sprintf("* Extracting E-OBS climatology for %s sites\n ", length(sites)))

  # download or read data
  eobs_data = lapply(c("tg","rr","elev"),function(x){
    # filename
    filename = sprintf("%s_%sdeg_reg_v15.0.nc",
                       x,
                       resolution)

    # if the file exist use the local file
    if (file.exists(sprintf("%s/%s", eobs_path, filename))){
      r = raster::brick(sprintf("%s/%s", eobs_path, filename))
      return(r)
    } else {
      stop('No E-OBS files found in the referred path !')
    }
  })

  # extract the yday and year strings and convert to numbers
  yday = as.numeric(format(as.Date(eobs_data[[1]]@z$Date),"%j"))
  years = as.numeric(format(as.Date(eobs_data[[1]]@z$Date),"%Y"))

  # construct validation data using the helper function
  # format_data() above
  validation_data = lapply(sites, function(x) {
    cat(sprintf(" |_ processing: %s of %s sites\r",
                which(sites == x),
                length(sites)))
    format_data(site = x,
                offset = offset)
  })

  # add proper list names (this should be the same name as the
  # site name added by the format_data() helper function)
  names(validation_data) = sites

  # screen for NULL values (climate grid cells with NA - near water bodies)
  validation_data = validation_data[!unlist(lapply(validation_data,is.null))]

  # return the formatted data
  return(validation_data)
}

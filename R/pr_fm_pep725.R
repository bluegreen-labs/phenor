#' Preprocessing of PEP725 data into a format which can be ingested
#' by the optimization routines.
#'
#' Some pre-processing steps are required as downloading the PEP725.
#' So, for the species of interest download the separate zipped files.
#' Unzip all files and put the respective scientific data in one folder.
#'
#' The routine requires E-OBS climate, which can be downloaded from
#' the E-OBS website and should be put unzipped in a single folder.
#' (http://www.ecad.eu/download/ensembles/ensembles.php).
#'
#' @param pep_path path to the PEP725 data (species files only)
#' @param eobs_path path to regular grid E-OBS data.
#' @param bbch which phenophase (bbch) to use (default = 11)
#' @param species species to select from merged PEP725 file
#' @param offset offset of the time series in DOY (default = 264, sept 21)
#' @param count minimum number of acquisitions per location
#' @param resolution resolution of the E-OBS data (0.25 or 0.5, default = 0.25)
#' @param pep_data prefiltered subset of the merged PEP725 data obtained through
#' merge_PEP725(), else the complete data contained in pep_path will be used
#' @return returns a nested list of site locations, their respective
#' phenological metrics and matching environmental data as extracted from
#' the E-OBS product (corrected for altitude using a lapse rate of 5C/km.)
#' @keywords phenology, model, preprocessing
#' @export
#' @examples
#'
#' # run with default settings
#' # looks for transition date files derived
#' # through phenocamr in your home directory
#' # change the path to match your setup
#' \dontrun{
#' phenocam_data = format_pep725(pep_path = "~/pep_data/",
#'                               eobs_path = "~/eobs_data/")
#'}

pr_fm_pep725 <- function(
  pep_path = tempdir(),
  eobs_path = tempdir(),
  bbch = "11",
  species = NULL,
  offset = 264,
  count = 60,
  resolution = 0.25,
  pep_data
  ) {

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
    tmean = raster::extract(eobs_data[[1]], points)
    tmin = raster::extract(eobs_data[[4]], points)
    tmax = raster::extract(eobs_data[[5]], points)

    # close to water bodies raster cells can be empty and return
    # NA, in this case set the output to NULL and remove these
    # values later on
    if (all(is.na(tmean))){
      return(NULL)
    } else {
      precipitation = raster::extract(eobs_data[[2]], points)
      tmean = matrix(rep(tmean, nrow(pep_subset)),
                           length(tmean),
                           nrow(pep_subset))
      tmin = matrix(rep(tmin, nrow(pep_subset)),
                           length(tmean),
                           nrow(pep_subset))
      tmax = matrix(rep(tmax, nrow(pep_subset)),
                           length(tmean),
                           nrow(pep_subset))
      precipitation = matrix(rep(precipitation, nrow(pep_subset)),
                             length(precipitation),
                             nrow(pep_subset))
    }

    # calculate long term mean temperature
    ltm = as.vector(unlist(by(tmean[which(years >= 1980),1],
             INDICES = yday[which(years >= 1980)],
             mean,
             na.rm = TRUE)))[1:365]

    # calculate lapse rate and correct temperatures
    lapse_rate = as.vector(
      (raster::extract(eobs_data[[3]],
                       points[1]) - pep_subset$alt[1]) * 0.005)
    tmean = rbind(tmean + lapse_rate, pep_subset$year)
    tmin = rbind(tmin + lapse_rate, pep_subset$year)
    tmax = rbind(tmax + lapse_rate, pep_subset$year)
    ltm = ltm + lapse_rate

    # now for all years create subsets for temperature
    # and precipitation (wrap this in a function and a do call)
    Ti = apply(tmean, 2, function(x){
      layers = which((years == (x[length(x)] - 1) & yday >= offset) |
                       (years == x[length(x)] & yday < offset))[1:365]
      return(x[layers])
    })

    Tmini = apply(tmin, 2, function(x){
      layers = which((years == (x[length(x)] - 1) & yday >= offset) |
                       (years == x[length(x)] & yday < offset))[1:365]
      return(x[layers])
    })

    Tmaxi = apply(tmax, 2, function(x){
      layers = which((years == (x[length(x)] - 1) & yday >= offset) |
                       (years == x[length(x)] & yday < offset))[1:365]
      return(x[layers])
    })

    Pi = apply(rbind(precipitation, pep_subset$year), 2, function(x){
      layers = which((years == (x[length(x)] - 1) & yday >= offset) |
                       (years == x[length(x)] & yday < offset))[1:365]
      return(x[layers])
    })

    # shift data when offset is < 365
    if (offset < 365){
      doy_neg = c((offset - 366):-1,1:(offset - 1))
      doy = c(offset:365,1:(offset - 1))
    } else {
      doy = doy_neg = 1:365
    }

    # calculate daylength using the doy
    l = ncol(Ti)
    Li = daylength(doy = doy, latitude = pep_subset$lat[1])
    Li = matrix(rep(Li,l),length(Li),l)

    # recreate the validation data structure (new format)
    # but with concatted data
    data = list("site" = site,
                "location" = c(pep_subset$lat[1], pep_subset$lon[1]),
                "doy" = doy_neg,
                "ltm" = ltm,
                "transition_dates" = pep_subset$day,
                "year" = pep_subset$year,
                "Ti" = Ti,
                "Tmini" = Tmini,
                "Tmaxi" = Tmaxi,
                "Li" = Li,
                "Pi" = Pi,
                "VPDi" = NULL,
                "georeferencing" = NULL
    )

    # return the list
    return(data)
  }

  cat("* Merging and cleaning PEP725 data files in: \n")
  cat(sprintf("  %s\n", pep_path))

  # User may provide prefiltered merge_PEP725 dataset
  if (missing(pep_data)){
    pep_data = merge_pep725(path = pep_path)
    }

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

  cat(sprintf(" |_ excluding sites with: < %s site years of observations \n", count))
  # count number of observations for each
  # unique pep725 site, this does not screen
  # for them being continuous or not !
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

  # read E-OBS data
  eobs_data = lapply(c("tg","rr","elev","tn","tx"),function(x){
    # filename
    filename=list.files(eobs_path,sprintf("%s_%sdeg_reg[^/]*\\.nc",x,resolution))
    # if the file exist use the local file
    if (length(filename)>0){
      r = raster::brick(sprintf("%s/%s", eobs_path, filename))
      return(r)
    } else {
      stop('No E-OBS files found in the referred path !')
    }
  })

  # extract the yday and year strings and convert to numbers
  yday = as.numeric(format(as.Date(eobs_data[[1]]@z$Date),"%j"))
  years = as.numeric(format(as.Date(eobs_data[[1]]@z$Date),"%Y"))

  # track progress
  pb = utils::txtProgressBar(min = 0, max = length(sites), style = 3)
  env = environment()
  i = 0

  # process data
  validation_data = lapply(sites, function(x) {
    utils::setTxtProgressBar(pb, i + 1)
    assign("i", i+1, envir = env)
    format_data(site = x,
                offset = offset)
  })

  # close progress bar
  close(pb)

  # add proper list names (this should be the same name as the
  # site name added by the format_data() helper function)
  names(validation_data) = sites

  # assign a class for post-processing
  class(validation_data) = "phenor_time_series_data"

  # screen for NULL values (climate grid cells with NA - near water bodies)
  validation_data = validation_data[!unlist(lapply(validation_data,is.null))]

  # return the formatted data
  return(validation_data)
}

#' phenocam_DB
#'
#' Data structure containing budburst dates for all deciduous broadleaf
#' (DB) sites included in the first PhenoCam data release. These data
#' are described in detail in Richardson et al. (2017)
#'
#' @format A structured list with 6 variables:
#' \describe{
#'   \item{site}{PhenoCam site name}
#'   \item{location}{geographic location}
#'   \item{doy}{DOY (1 - 365 as a vector)}
#'   \item{transition_dates}{budburst dates, as per 25 percent of the Gcc amplitude }
#'   \item{Ti}{mean daily temperature, matrix 365 rows x # site years columns}
#'   \item{Li}{daylength, matrix 365 rows x # site years columns}
#' }
"phenocam_DB"

#' phenocam_EN
#'
#' Data structure containing budburst dates for all evergreen needleleaf
#' (EN) sites included in the first PhenoCam data release. These data
#' are described in detail in Richardson et al. (2017)
#'
#' @format A structured list with 6 variables:
#' \describe{
#'   \item{site}{PhenoCam site name}
#'   \item{location}{geographic location}
#'   \item{doy}{DOY (1 - 365 as a vector)}
#'   \item{transition_dates}{budburst dates, as per 25 percent of the Gcc amplitude }
#'   \item{Ti}{mean daily temperature, matrix 365 rows x # site years columns}
#'   \item{Li}{daylength, matrix 365 rows x # site years columns}
#' }
"phenocam_EN"

#' phenocam_GR
#'
#' Data structure containing budburst dates for all grassland (GR)
#' sites included in the first PhenoCam data release. These data
#' are described in detail in Richardson et al. (2017)
#'
#' @format A structured list with 6 variables:
#' \describe{
#'   \item{site}{PhenoCam site name}
#'   \item{location}{geographic location}
#'   \item{doy}{DOY (1 - 365 as a vector)}
#'   \item{transition_dates}{budburst dates, as per 25 percent of the Gcc amplitude }
#'   \item{Ti}{mean daily temperature, matrix 365 rows x # site years columns}
#'   \item{Li}{daylength, matrix 365 rows x # site years columns}
#' }
"phenocam_GR"

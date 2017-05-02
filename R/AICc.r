#' Calculates AICc values for a set of measured and predicted values
#' together with the number of model parameters used
#'
#' @param measured: a vector with measurement values to smooth
#' @param predicted: a vector with dates / time steps
#' @param k: optional values to weigh the loess fit with
#' @keywords model selection, AIC, Akaike's Information Criterion
#' @export
#' @examples
#'
#' \dontrun{
#'
#' model_AIC = AICc(measured, predicted, k)
#'
#' }

# custom AIC function which accepts loess regressions
AICc = function(measured, predicted, k){

  # calculate number of observations
  n = length(measured)

  # calculatue residual sum of squares
  RSS = sum((measured - predicted)^2)

  # AIC
  AIC = 2*k + n * log(RSS/n)

  # AICc
  AICc  = AIC + (2 * k * (k + 1)) / (n - k - 1)

  # return both AIC
  return(list("AIC"=AIC,"AICc"=AICc))
}

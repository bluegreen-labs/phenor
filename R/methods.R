#' Plotting method for fit model output
#'
#' @param data input data generated using the pr_fit() function
#' @return a plot with fit statistics
#' @keywords phenology, model, accuracy
#' @export

# plot data if requested
plot.pr_fit <- function(data){
  stopifnot(inherits(data, "pr_fit"))
  plot(data$measured,
       data$predicted,
       main = data$model,
       xlab = "onset DOY Measured",
       ylab = "onset DOY Modelled",
       pch = 19,
       tck = 0.02)
  graphics::abline(0,1)
  graphics::legend("topleft",
                   legend = sprintf("RMSE: %s",
                                    round(data$rmse,2)),bty='n')
  graphics::legend("top",
                   legend = sprintf("RMSE NULL: %s",
                                    round(data$rmse_null,2)),bty='n')
  graphics::legend("bottomright",
                   legend = sprintf("AICc: %s",
                                    round(data$aic$AICc)),bty='n')
}

#' Print summary values for fit model output
#'
#' @param data input data generated using the pr_fit() function
#' @return a plot with fit statistics
#' @keywords phenology, model, accuracy
#' @export

summary.pr_fit <- function(data){
  stopifnot(inherits(data, "pr_fit"))

  # read parameter ranges
  d <- t(pr_parameters(model = data$model))
  d <- cbind(rownames(d),d, as.vector(round(data$par,2)))
  colnames(d) <- c("name","lower", "upper", "best")

  # basic statistics
  title <- paste0("Statistics for the ", data$mode, " model")
  message(title)
  message(paste(c(rep("-", nchar(title)),"\n"), collapse = ""))

  message(paste0(" - Model accuracy (RMSE): ", round(data$rmse, 2)))
  message(paste0(" - NULL model (RMSE): ",
                 round(data$rmse_null, 2)))
  message(paste0(" - AIC value : ",
                 round(data$aic[[1]], 2),"\n"))

  message("Parameter values")
  message(paste(c(rep("-", nchar(title)),"\n"), collapse = ""))
  message(paste("", rownames(d), collapse = "\t"))
  message(paste("", round(data$par, 2), collapse = "\t"))
  return(d)
}




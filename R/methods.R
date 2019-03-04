# plot data if requested
plot.pr_fit <- function(data){
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
                                    round(data$rmse)),bty='n')
  graphics::legend("top",
                   legend = sprintf("RMSE NULL: %s",
                                    round(data$rmse_null)),bty='n')
  graphics::legend("bottomright",
                   legend = sprintf("AICc: %s",
                                    round(data$aic$AICc)),bty='n')
}

# print fit data
print.pr_fit <- function(data){
  # basic statistics
  print(data$rmse)
  print(data$rmse_null)
  print(data$aic)

  # print summary statistics
  print(summary(stats::lm(data$measured ~ data$predicted)))
}




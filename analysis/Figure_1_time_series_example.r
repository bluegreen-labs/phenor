# Figure 1.
library(phenor)

# read in the demo time series contained in the package,
# alternatively the data can be found in the Nature Scientific Data
# publication cited on the phenor github page
time_series = read.table(file.path(path.package("phenor"),
                                   "extdata/bartlett_DB_0003_3day.csv"),
                         header = TRUE,
                         sep = ",")

transition_dates = read.table(file.path(path.package("phenor"),
                                        "extdata/bartlett_DB_0003_3day_transition_dates.csv"),
                         header = TRUE,
                         sep = ",")

# subset phenophase data per "season"
# could probably be done with ggplot magrittr pipe magic, but hey cut
# me some slack
rising = transition_dates[transition_dates$gcc_value=="gcc_90" &
                            transition_dates$direction == "rising",]
falling = transition_dates[transition_dates$gcc_value=="gcc_90" &
                             transition_dates$direction == "falling",]

# set device
pdf("~/Figure_1_time_series_example.pdf", 14, 4)
par(tck = 0.02)
plot(
  as.Date(time_series$date),
  time_series$gcc_90,
  type = 'n',
  xlab = "Year",
  ylab = "Gcc"
)
lines(
  as.Date(time_series$date),
  time_series$smooth_gcc_90,
  col = "grey",
  lwd = 3
)
points(as.Date(time_series$date),
       time_series$gcc_90,
       pch = 19)

# plot transition dates
points(
  as.Date(rising$transition_25),
  rising$threshold_25,
  col = "blue",
  pch = 15,
  cex = 2
)
points(
  as.Date(falling$transition_25),
  falling$threshold_25,
  col = "red",
  pch = 19,
  cex = 2
)
dev.off()

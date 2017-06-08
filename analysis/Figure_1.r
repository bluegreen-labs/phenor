# Figure 1.

# An example Gcc time series.
# this code should be run within the directory which contains the code
# set your working directory accordingly.

# read in the model data structure (contains the geographic location)
time_series = read.table("../data/data_record_5/bartlett_DB_0003_3day.csv",
                         header = TRUE,
                         sep = ",")

transition_dates = read.table("../data/data_record_6/bartlett_DB_0003_3day_transition_dates.csv",
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
pdf("../output/Figure_1.pdf", 12, 5)
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
  rising$threshold_50,
  col = "green",
  pch = 15,
  cex = 2
)
points(
  as.Date(falling$transition_25),
  falling$threshold_50,
  col = "brown",
  pch = 19,
  cex = 2
)
dev.off()

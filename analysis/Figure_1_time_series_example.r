# Figure 1.

# make sure to install the required libraries
# and load them
library(phenor)
library(ggplot2)
library(ggthemes)

# list required files
data_dir = sprintf("%s/extdata/figure_1/",file.path(path.package("phenor")))
time_series_files = list.files(data_dir,
                         "*3day.csv",
                         full.names = TRUE)
transition_dates = list.files(data_dir,
                              "*transition_dates.csv",
                              full.names = TRUE)

# merge data files
time_series = do.call("rbind",
                      lapply(time_series_files,function(x){
  df = read.table(x,header = TRUE, sep = ",")
  df$date = as.Date(df$date)
  start = df$date[1]
  end = df$date[nrow(df)]
  df = df[df$date >= start+90 & df$date <= end-90,]
  filename = unlist(strsplit(basename(x),"_"))
  df$site = filename[1]
  df$veg_type = filename[2]
  df$roi = filename[3]
  return(df)
}))

transition_dates = do.call("rbind",lapply(transition_dates,function(x){
  read.table(x,header = TRUE, sep = ",")
}))

# creaete subset based upon desired time series
transition_dates = transition_dates[transition_dates$gcc_value == "gcc_90",]

# merge time_series with transition_data files

# create empty slots in time_series file
time_series$transition_25 = NA
time_series$threshold_25 = NA
time_series$direction = NA

# create a temporary file with similar dimensions
tmp = as.data.frame(matrix(NA,nrow(transition_dates),ncol(time_series)))
colnames(tmp) = colnames(time_series)
tmp$direction = transition_dates$direction

# copy required transition date data
tmp$transition_25 = transition_dates$transition_25
tmp$threshold_25 = transition_dates$threshold_25
tmp$site = transition_dates$sitename

# bind the data together
combined_data = rbind(time_series,tmp)

# convert dates to proper format
combined_data$transition_25 = as.Date(combined_data$transition_25)
combined_data$date = as.Date(combined_data$date)

# ggplot plotting routines
p = ggplot(data = combined_data) +
  geom_line(aes(x = date, y = smooth_gcc_90),
            lwd = 1.2) +
  geom_point(aes(x = transition_25,
                 y = threshold_25,
                 colour = direction,
                 shape = direction),
             size = 2) +
  facet_wrap(~ site, ncol = 2) +
  theme_hc() +
  theme(strip.text.x = element_text(size = 12),
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=12, vjust = 2),
        axis.title.y = element_text(size=12),
        legend.position="none") +
  scale_color_manual(values = c("red", "blue")) +
  scale_shape_manual(values = c(19,15)) +
  xlab("Date") +
  ylab("Gcc")

# export figure (to home directory)
pdf("~/Figure_1_time_series_example.pdf", 14, 10)
  plot(p)
dev.off()

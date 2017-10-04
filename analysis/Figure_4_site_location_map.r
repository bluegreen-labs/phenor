# Figure 4.

# An overview map of the PhenoCam deciduous broadleaf sites as included
# in the PhenoCam Data Paper (Richardson et al. 2017).
# this code should be run within the directory which contains the code
# set your working directory accordingly.

# load libraries
library(ggplot2)
library(maptools)
library(phenor)

# read in the model data structure (contains the geographic location)
data("phenocam_DB")
data("phenocam_EN")
data("phenocam_GR")

# flatten the format
phenocam_DB = flat_format(phenocam_DB)
phenocam_EN = flat_format(phenocam_EN)
phenocam_GR = flat_format(phenocam_GR)

# create a spatial points object and transform to polar coords
sites_DB = unique(data.frame(t(phenocam_DB$location[1:2,])))
sites_EN = unique(data.frame(t(phenocam_EN$location[1:2,])))
sites_GR = unique(data.frame(t(phenocam_GR$location[1:2,])))

colnames(sites_GR) = c("lat","lon")
colnames(sites_EN) = c("lat","lon")
colnames(sites_DB) = c("lat","lon")

p = ggplot() +
  geom_polygon(data = map_data("world"),
               aes(x=long, y=lat, group=group),
               fill='grey20',
               colour = 'white',
               size = 0.1) +
  scale_x_continuous(breaks=seq(-160,-50,20)) +
  scale_y_continuous(breaks=seq(23,65,10)) +
  coord_cartesian(xlim = c(-160, -50), ylim = c(23, 65)) +
  geom_point(data = sites_DB,
             aes(x = lon, y = lat),
             col = "yellow",
             alpha = 0.8,
             pch = 15,
             size = 3) +
  geom_point(data = sites_GR,
             aes(x = lon, y = lat),
             col = "lightblue",
             alpha = 1,
             pch = 17,
             size = 3) +
  geom_point(data = sites_EN,
             aes(x = lon, y = lat),
             col = "red",
             alpha = 1,
             pch = 19,
             size = 3) +
  xlab(expression(paste("Longitude (", degree, ")"))) +
  ylab(expression(paste("Latitude (", degree, ")"))) +
  theme_bw()

# export device
pdf("~/Figure_2_site_location_map.pdf",7.5,4.5)
  plot(p)
# device off
dev.off()

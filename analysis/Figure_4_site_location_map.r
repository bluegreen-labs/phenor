# Figure 4.

# An overview map of the PhenoCam deciduous broadleaf sites as included
# in the PhenoCam Data Paper (Richardson et al. 2017).
# this code should be run within the directory which contains the code
# set your working directory accordingly.

# load libraries
library(maps)
library(mapdata)
library(maptools)
#library(sf)
library(phenor)

setwd("/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/")

# set coordinate systems to use
lat_lon =  CRS("+init=epsg:4326") # lat lon
us =  CRS("+init=epsg:2163") # us projection

# read in the model data structure (contains the geographic location)
data("phenocam_DB")
data("phenocam_EN")
data("phenocam_GR")

# flatten the format
phenocam_DB = flat_format(phenocam_DB)
phenocam_EN = flat_format(phenocam_EN)
phenocam_GR = flat_format(phenocam_GR)

# download the blue marble data if it doesn't
# exist
if (!file.exists("blue_marble.tif")) {
  download.file("http://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=526312&cs=rgb&format=TIFF&width=5400&height=2700",
                "blue_marble.tif")
}

# read in the raster map and
# set the extent, crop to extent and reproject to polar
r = raster::brick("blue_marble.tif")
e = raster::extent(c(-160,-50,23,65))
r_crop = raster::crop(r,e)

# traps NA values and sets them to 1
r_crop[is.na(r_crop)] = 1

# create a spatial points object and transform to polar coords
sites_DB = SpatialPoints(t(phenocam_DB$location[2:1,]),lat_lon)
sites_EN = SpatialPoints(t(phenocam_EN$location[2:1,]),lat_lon)
sites_GR = SpatialPoints(t(phenocam_GR$location[2:1,]),lat_lon)

# grab country outlines
countries = map('worldHires',
                c('Canada','USA','Mexico'),
                fill=TRUE,
                plot = FALSE)

# grab IDs
IDs = sapply(strsplit(countries$names, ":"), function(x) x[1])

# convert to a spatial object, in lat lon
countries = map2SpatialLines(countries,
                             IDs = IDs,
                             proj4string = lat_lon)

# export device
pdf("../output/Figure_2.pdf",14.5,8.5)

# margin settings etc
par(
  oma = c(5, 6.5, 5, 6.5),
  bg = "white",
  fg = "white",
  tck = 0.01,
  lwd = 2,
  cex.axis = 1,
  cex = 1.3,
  col.axis = "grey25",
  col.lab = "grey25",
  col.main = "grey25",
  col.sub = "grey25",
  col = "grey25"
)

# plot the blue marble raster
raster::plotRGB(r_crop)

# plot grid lines
grid(lwd=1,
     lty=2)

# plot continent outlines
plot(countries,
     col = "grey50",
     lwd = 0.8,
     add = TRUE)

# plot sites
points(sites_DB,
       col = "yellow",
       pch = 0,
       cex = 1.2)

points(sites_EN,
       col = "red",
       pch = 1,
       cex = 1.2)

points(sites_GR,
       col = "lightblue",
       pch = 2,
       cex = 1.2)

# plot axis and axis labels
axis(
  1,
  col.ticks = "grey25"
)
axis(
  2,
  col.ticks = "grey25"
)

# labels
mtext(expression(paste("Longitude (", degree, ")")),
      1,
      3,
      cex = 1.3)
mtext(expression(paste("Latitude (", degree, ")")),
      2,
      3,
      cex = 1.3)

# add a fat box around the map
box(lwd = 2)

# device off
dev.off()

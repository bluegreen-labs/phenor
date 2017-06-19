# Create spatial representations from scratch
# for full transparency in processing. Take
# into account that this might take a while
# for regenerating the maps as shown in
# the original manuscript. Mainly downloading
# and processing tiled Daymet data is rather slow.

# load the necessary libraries
library(phenor)
library(daymetr)
library(maps)
library(RColorBrewer)

# path to store phenor data
path = "~/phenor_files"

# check if the path exists if not create it
if ( !dir.exists(path) ){
  dir.create(path)
}

# set working directory, create path parameters in
# necessary functions (TODO)
setwd(path)

# # download daymet tiles and format the data
# # this might take a while (be warned), first
# # list all the tiles needed
# tiles = c(
#   11936, 12295, 12296, 12297,
#   12114, 12115, 12116, 12117,
#   11934, 11935, 11754, 11755,
#   11756, 12294, 11937)
#
# # download all daymet tiles, include the preceeding year
# # as we need this data too (see start_yr).
# daymetr::download_daymet_tiles(tiles = tiles,
#                                param = c("tmin","tmax"),
#                                start_yr = 2010,
#                                end_yr = 2011)
#
# # now calculate the mean daily temperature from tmin and tmax
# lapply(tiles, function(x)daymetr::daymet_tmean(path = path,
#                                                tile = x,
#                                                year = 2010,
#                                                internal = FALSE))
#
# lapply(tiles, function(x)daymetr::daymet_tmean(path = path,
#                                                tile = x,
#                                                year = 2011,
#                                                internal = FALSE))
#
# # now format all tiles acoording to the phenor format
# # (default settings)
# lapply(tiles, function(x)daymet_data = format_daymet_tiles(
#   path = path,
#   tile = x,
#   year = 2011,
#   internal = FALSE))
#
# # get cmip5 data for the end of the century
cmip5_data_2100 = format_cmip5(year = 2100)
cmip5_data_2011 = format_cmip5(year = 2011)

# download berkeley earth data (global data but subset for CONUS)
# for the same year as the Daymet data (scale comparison)
#be_data = format_berkeley_earth(year = 2011)

# create a set of demo parameters
# optimize for default Deciduous Broadleaf
# data (included in the phenor package)
# I use the simplest Thermal Time model (TT)
par = model_validation(model = "TT",
                       control = list(max.call = 5000,
                                      temperature = 10000))

# create maps using estimated parameters and
# the estimate_phenology routine
cmip5_map_2100 = estimate_phenology(par = par,
                               model = "TT",
                               data = cmip5_data_2100)

# create maps using estimated parameters and
# the estimate_phenology routine
cmip5_map_2011 = estimate_phenology(par = par,
                                    model = "TT",
                                    data = cmip5_data_2011)

#be_map = estimate_phenology(par = par,
#                            model = "TT",
#                            data = be_data)

daymet_map = estimate_phenology(par = par,
                                model = "TT",
                                path = path)

# Generate the final plot comparing model output
# across various scales and time frames use
# a pdf device for quality graphs
pdf("~/Figure_5.pdf",12,14)

# set margins and general layout
# of subplots
par(oma = c(5,4,5,2))
layout(matrix(c(1,1,2,2,
                1,1,2,2,
                4,4,2,2,
                4,4,2,2,
                3,3,5,5,
                3,3,6,6),
              6, 4, byrow = TRUE))

# define the colours to use
cols = colorRampPalette(brewer.pal(9,'RdBu'))(100)
zlim =  c(90,190)

# what data do we want to map (generated above)
mapdata = c("cmip5_map_2011","daymet_map","cmip5_map_2100")

# loop over all datasets and plot them in the
# correct order with adjusted settings
for (i in 1:3){
  if (i != 3){
    par(mar = c(1,1,1,1),
        cex = 1.1,
        cex.lab = 1.1)
  } else {
    par(mar = c(0,1,2,1))
  }

  image(get(mapdata[i]),
                xlab = '',
                ylab = '',
                xaxt = 'n',
                yaxt  = 'n',
                tck = 0.02,
                bty = 'n',
                zlim = zlim,
                col = cols)
  box()
  map("world",
      col = "grey25",
      add = TRUE)

  #  Headers to the subplots
  if (i == 1){
    e = as.vector(extent(daymet_map))
    rect(e[1],e[3],e[2]-0.4,e[4],
         lty = 1,
         lwd = 1.5,
         border = "black")
    mtext("CMIP5 (IPSL-CM5A-MR, 2011)",
          side = 3,
          line = 1,
          cex = 1.3)
    axis(1,
         tck = 0.02,
         labels = FALSE)
    axis(2,
         tck = 0.02)
    legend("bottomright",
           "a",
           bty = "n",
           cex = 1.5)

    # overplot the first plot
    # to insert a pretty arrow
    # to the plot to the right
    par(new = TRUE,
        mar = c(0,0,0,0))
    plot(0,
         type = 'n',
         xlim = c(0,1),
         ylim = c(0,1),
         axes = FALSE,
         bty = 'n')
    segments(0.95,
             0.55,
             0.95,
             0.4,
             lwd = 1.5)
    arrows(0.95,0.4,1.04,0.4,
           angle = 30,
           length = 0.1,
           lwd = 1.5)

  } else if ( i == 3 ){
    mtext("CMIP5 (IPSL-CM5A-MR, 2100)",
          side = 3,
          line = 1,
          cex = 1.3)
    axis(1,
         tck = 0.02)
    axis(2,
         tck = 0.02)
    mtext(text = "longitude",
          side = 1,
          line = 3,
          cex = 1.5)
    legend("bottomright",
           "d",
           bty = "n",
           cex = 1.5)
  } else {
    mtext("Daymet (2011)",
          side = 3,
          line = 1,
          cex = 1.3)
    axis(1,
         tck = 0.02)
    axis(4,
         tck = 0.02)
    mtext(text = "longitude",
          side = 1,
          line = 3,
          cex = 1.5)
    legend("bottomright",
           "c",
           bty = "n",
           cex = 1.5)

    # some locations annotated for reference
    # showing the fine grained structure in the
    # spatial Daymet data
    segments(-68.921274, 45.904354, -68, 44,
             lwd = 1,
             lty = 2)
    text(-68, 43.8, "Mt. Katahdin", cex = 1.3)
    points(-71.063611,42.358056,
           pch = 19,
           cex = 1.5)
    text(-71.063611,
         42.358056,
         "Boston",
         cex = 1.3, pos = 2)


  }
}

# select those land cover pixels
# with more than 1/4 of the pixel covered
# igbp classes 1/4/5/10 are included in the
# package as igbp_#
m = igbp_4 > 0.25

# plot a filtered version of the be_map data
par(mar = c(1,1,1,1))
image(cmip5_map_2011 * m,
      xlab = '',
      ylab = '',
      xaxt = 'n',
      yaxt  = 'n',
      tck = 0.02,
      bty = 'n',
      zlim = zlim,
      col = cols)
box()
map("world",
    col = "grey25",
    add = TRUE)
axis(1,
     tck = 0.02,
     labels = FALSE)
axis(2,
     tck = 0.02)
mtext(text = "latitude",
      side = 2,
      line = 3,
      cex = 1.5)
legend("bottomright",
       "b",
       bty = "n",
       cex = 1.5)

# colour legend
par(mar = c(3,1,5,1))
imageScale(cmip5_map_2011,
           col = cols,
           zlim = zlim,
           cex = 1.5,
           axis.pos = 1)

mtext("DOY",
      1,
      3,
      cex = 1.5)

# close device to generate a valid output file
dev.off()

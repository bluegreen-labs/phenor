# Figure 3

# showing the different PFT datasets and model formulations
# and their respective model skill

# load the data
comparison = readRDS("/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/comparison.rda")
melaas = readRDS("/data/Dropbox/Research_Projects/working/phenocam_model_comparison/data/melaas.rda")

# set the outline
par(mfrow=c(4,1), oma = c(5,0,0,0))

# plot all panels
par(mar=c(0,5,2,1))
plot_model_comparison(comparison$phenocam_DB, names = FALSE)
legend("topleft", "Deciduous Broadleaf", bty = "n", cex = 1.5)
legend("topright", "a)", bty = "n", cex = 2)

par(mar=c(0,5,1,1))
plot_model_comparison(melaas, names = FALSE)
legend("topleft", "Deciduous Broadleaf (Melaas et al. 2016)", bty = "n", cex = 1.5)
legend("topright", "b)", bty = "n", cex = 2)

par(mar=c(1,5,1,1))
plot_model_comparison(comparison$phenocam_EN, names = FALSE)
legend("topleft", "Evergreen Needleleaf", bty = "n", cex = 1.5)
legend("topright", "c)", bty = "n", cex = 2)

par(mar=c(1,5,0,1))
plot_model_comparison(comparison$phenocam_GR, names = TRUE)
legend("topleft", "Grassland", bty = "n", cex = 1.5)
legend("topright", "d)", bty = "n", cex = 2)

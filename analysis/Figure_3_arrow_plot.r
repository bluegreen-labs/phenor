# Figure 3.
#
# Arrow plot comparing two model optimizations
# and the difference in model output.
# Mean differences are used across different runs within
# a given model.
library(phenor)

# quick comparison with default settings (1 random seed)
comparison = model_comparison(models = c("TT","PTT"),
                              random_seeds = 1,
                              control = list(max.call = 40000,
                                             temperature = 10000))
# plot the data nicely
pdf("~/Figure_3_arrow_plot.pdf",7,5)
arrow_plot(comparison)
dev.off()

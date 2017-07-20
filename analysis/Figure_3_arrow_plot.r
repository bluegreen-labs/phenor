# Figure 3.
#
# Arrow plot comparing two model optimizations
# and the difference in model output.
# Mean differences are used across different runs within
# a given model.
library(phenor)

# quick comparison with default settings (1 random seed)
pdf("~/Figure_3_arrow_plot.pdf",7,5)
comparison = model_comparison(models = c("TT","PTT"),
                              random_seeds = 1,
                              control = list(max.call = 10000,
                                             temperature = 10000))
arrow_plot(comparison)
dev.off()

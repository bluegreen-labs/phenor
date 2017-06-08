# Figure 4.
#
# Arrow plot comparing two model optimizations
# and the difference in model output.
# Mean differences are used across different runs within
# a given model.
library(phenor)

# quick comparison with default settings (3 random seeds)
comparison = model_comparison(models = c("TT","PTT"))
arrow_plot(comparison)

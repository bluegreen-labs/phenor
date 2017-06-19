# Figure 3.
#
# Model validation plot as shown when using
# the model validion function.
# Default parameters are used

library(phenor)

pdf("~/Figure_3.pdf",7,5)
model_validation(control = list(max.call = 10000,
                                temperature = 10000))
dev.off()

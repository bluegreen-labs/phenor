# Figure 2.
#
# Model validation plot as shown when using
# the model validion function.
# Default parameters are used

library(phenor)

pdf("~/Figure_2_model_validation.pdf",7,5)
model_validation(control = list(max.call = 10000,
                                temperature = 10000))
dev.off()

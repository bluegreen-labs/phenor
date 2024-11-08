try(detach("package:phenor", unload = TRUE))
library(phenor)

set.seed(0)

par_bt = phenor::pr_fit(
  data = phenocam_DB,
  model = "TT",
  method = "bayesiantools",
  control = list(
    sampler = "DEzs",
    settings = list(
      burnin = 10,
      iterations = 1000
      )
  )
)

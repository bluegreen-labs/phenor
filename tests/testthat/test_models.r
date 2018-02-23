# Phenor unit tests

# test all phenology models
test_that("test model runs",{

  # load parameter ranges
  models = utils::read.table(sprintf("%s/extdata/parameter_ranges.csv",
                                         path.package("phenor")),
                                 header = TRUE,
                                 sep = ",",
                                 stringsAsFactors = FALSE)$model

  # test models
  model_runs = try(model_comparison(random_seeds = 1,
                                    models = models,
                                    control = list(max.call = 10)))

  # see if any of the runs failed
  check = !inherits(model_runs, "try-error")

  # check if no error occured
  expect_true(check)
})

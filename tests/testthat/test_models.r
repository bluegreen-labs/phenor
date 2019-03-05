# Phenor unit tests

# test all phenology models
test_that("test model runs",{

  # load parameter ranges
  models <- unique(
    utils::read.table(
      system.file(
        "extdata",
        "parameter_ranges.csv",
        package = "phenor",
        mustWork = TRUE
      ),
      header = TRUE,
      sep = ",",
      stringsAsFactors = FALSE
    )$model
  )

  # test models
  expect_output(pr_fit_compare(random_seeds = 1,
                                 models = models,
                                 control = list(max.call = 1)))
})

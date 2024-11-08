# Phenor unit tests

# test all phenology models
test_that("test models",{

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

  # subset models ignore new ones Alison
  models <- models[!grepl("W",models)]

  # test models
  expect_output(pr_fit_comparison(
    data = phenocam_DB[c(1:2)],
    random_seeds = 1,
    models = models[1:2],
    method = "bayesiantools",
    control = list(
      sampler = "DEzs",
      settings = list(
        burnin = 10,
        iterations = 1000)
      )
    )
  )

  # test models
  expect_output(pr_fit_comparison(
    data =  phenocam_DB[c(1:2)],
    random_seeds = 1,
    models = models[1:2],
    control = list(max.call = 5)
    )
  )

  # fit model
  fit <- pr_fit(control=list(max.call = 5))

  # fit model test
  expect_type(
    pr_fit(
      control=list(max.call = 5),
      cost = cvmae),
    "list"
  )

  expect_type(pr_calc_temp_sens(
      par = fit$par,
      data = phenocam_DB,
      model = "TT"
    ),
    "list"
  )

  # missing model
  expect_error(
    pr_calc_temp_sens(
      par = fit$par,
      data = phenocam_DB
      )
  )

  # missing input data
  expect_error(
    pr_calc_temp_sens(
      par = fit$par,
      model = "TT"
    )
  )

  # missing parameters
  expect_error(
    pr_calc_temp_sens(
      data = phenocam_DB,
      model = "TT"
    )
  )

  # return summary stats
  expect_message(summary(fit))
  expect_type(plot(fit),"list")

  # test CV
  expect_type(
    pr_cross_validate(
        k = 2,
        cv_seed = 123,
        models =  c("LIN"),
        data = phenocam_DB[c(1:2)]),
      "list"
    )
})

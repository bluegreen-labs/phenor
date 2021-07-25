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

  # test models
  expect_output(pr_fit_comparison(
    data =phenocam_DB[c(1:2)],
    random_seeds = 1,
    models = models,
    control = list(max.call = 5)))

  # test models
  expect_output(pr_fit_comparison(
    data =  phenocam_DB[c(1:2)],
    random_seeds = 1,
    models = models,
    method = "bayesiantools",
    control = list(
      sampler = "DEzs",
      settings = list(
        burnin = 10,
        iterations = 1000)
      )
    )
  )

  # fit model
  fit <- pr_fit(control=list(max.call = 5))
  fit_cvmae <- pr_fit(control=list(max.call = 5), cost = cvmae)
  expect_type(fit_cvmae, "list")

  temp_sens <- pr_calc_temp_sens(
    par = fit$par,
    data = phenocam_DB,
    model = "TT")

  expect_type(temp_sens, "list")

  expect_error(
    pr_calc_temp_sens(
      par = fit$par,
      data = phenocam_DB
      )
  )

  expect_error(
    pr_calc_temp_sens(
      par = fit$par,
      model = "TT"
    )
  )

  expect_error(
    pr_calc_temp_sens(
      data = phenocam_DB,
      model = "TT")
  )

  # return summary stats
  expect_message(summary(fit))
  expect_type(plot(fit),"list")

  # test CV
  l <- pr_cross_validate(
    k = 2,
    cv_seed = 123,
    models =  c("LIN"),
    data = phenocam_DB[c(1:2)])

  expect_type(l, "list")
})

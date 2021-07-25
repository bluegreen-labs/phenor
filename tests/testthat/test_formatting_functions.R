# test formatting functions

test_that("test formatting functions",{

  models <- c("LIN","PTT")

  # test models
  expect_output(pr_fit_comparison(random_seeds = 1,
                                  models = models,
                                  control = list(max.call = 10)))
})

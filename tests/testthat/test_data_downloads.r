# Phenor unit tests

# test all data downloads
test_that("test data downloads",{

  # download npn data
  expect_output(pr_dl_npn(species = 3,
                              path = tempdir(),
                              internal = FALSE))

  # download npn data internal
  expect_output(pr_dl_npn(species = 3,
                          internal = TRUE))

  # download berkeley earth data
  expect_output(pr_dl_be(year = 2011, path = tempdir()))
})

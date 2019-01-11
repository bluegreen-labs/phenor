# Phenor unit tests

# test all data downloads
test_that("test data downloads",{

  # download npn data
  expect_output(download_npn(species = 3,
                              path = tempdir(),
                              internal = FALSE))

  # download npn data internal
  expect_output(download_npn(species = 3,
                                       internal = TRUE))

  # download cmip5 data (TOO SLOW)
  # cmip_data = try (download_cmip5(year = 2011,
  #                                 path = tempdir(),
  #                                 model = "MIROC5",
  #                                 scenario = "rcp85"))

  # download berkeley earth data
  expect_output(download_berkeley_earth(year = 2011,
                                        path = tempdir()))
})

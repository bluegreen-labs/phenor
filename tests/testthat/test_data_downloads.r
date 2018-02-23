# Phenor unit tests

# test all data downloads
test_that("test data downloads",{

  # download npn data
  npn_data = try(download_npn(species = 3,
                              path = tempdir(),
                              internal = FALSE))

  # download npn data internal
  npn_data_internal = try(download_npn(species = 3,
                                       internal = TRUE))

  # download cmip5 data
  cmip_data = try (download_cmip5(year = 2011,
                                  path = tempdir(),
                                  model = "MIROC5",
                                  scenario = "rcp85"))

  # download berkeley earth data
  be_data = try(download_berkeley_earth(year = 2011,
                                        path = tempdir()))

  # see if any of the runs failed
  check = !inherits(npn_data, "try-error") &
          !inherits(npn_data_internal, "try-error") &
          !inherits(cmip_data, "try-error") &
          !inherits(be_data, "try-error")

  # check if no error occured
  expect_true(check)
})

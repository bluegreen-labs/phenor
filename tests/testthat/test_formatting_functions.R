# test formatting functions

# test all data downloads
test_that("test downloading functions",{

  # NPN data
  expect_message(
    pr_dl_npn(species = 3,
              internal = FALSE,
              start = "2001-01-01",
              end = "2001-12-31")
  )

  expect_error(
    pr_dl_npn(species = 3,
              internal = FALSE,
              start = "2000-01-01",
              end = "2000-12-31")
  )

  df <- pr_dl_npn(species = 3,
                  internal = TRUE,
                  start = "2001-01-01",
                  end = "2001-12-31")

  expect_type(df,"list")

})

test_that("test formatting functions",{

  df <- pr_dl_npn(species = 3,
                  internal = TRUE,
                  start = "2001-01-01",
                  end = "2001-12-31")

  df <- pr_fm_npn(df)

  expect_type(df,"list")

  # check npn meta-data function
  expect_output(check_npn_phenophases(list = TRUE))
  l <- check_npn_phenophases(phenophase = 61, list = FALSE)
  expect_type(l,"logical")
  l <- check_npn_species(species = 3, list = TRUE)
  expect_type(l,"list")
  l <- check_npn_species(list = TRUE)
  expect_type(l,"list")

})

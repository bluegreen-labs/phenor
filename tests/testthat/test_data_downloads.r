# Phenor unit tests

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

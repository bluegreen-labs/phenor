# test formatting functions

# test all data downloads
# test_that("test downloading functions",{
#
#   # NPN data
#   expect_message(
#     pr_dl_npn(species = 3,
#               internal = FALSE,
#               start = "2001-01-01",
#               end = "2001-12-31")
#   )
#
#   expect_error(
#     pr_dl_npn(species = 3,
#               internal = FALSE,
#               start = "2000-01-01",
#               end = "2000-12-31")
#   )
#
#   df <- pr_dl_npn(species = 3,
#                   internal = TRUE,
#                   start = "2001-01-01",
#                   end = "2001-12-31")
#
#   expect_type(df,"list")
#
# })

test_that("test formatting functions",{

  # df <- pr_dl_npn(
  #   species = 3,
  #   internal = TRUE,
  #   start = "2001-01-01",
  #   end = "2001-12-31"
  #   )
  #
  # df <- pr_fm_npn(df)
  #
  # expect_type(df,"list")

  # check npn meta-data function
  # expect_output(check_npn_phenophases(list = TRUE))
  # l <- check_npn_phenophases(phenophase = 61, list = FALSE)
  # expect_type(l,"logical")
  # expect_output(check_npn_phenophases(phenophase = 61, list = TRUE))
  # l <- check_npn_species(species = 3, list = TRUE)
  # expect_type(l,"list")
  # l <- check_npn_species(list = TRUE)
  # expect_type(l,"list")

  # csv data
  f <- data.frame(
    id = "site1",
    lat = 42,
    lon = -110,
    phenophase = "spring",
    year = 2000,
    doy = 120
    )

  t <- tempfile()
  write.table(
    f,
    t,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = ","
    )

  # format CSV data
  csv_data <- pr_fm_csv(t, phenophase = "spring")
  expect_type(csv_data, "list")

  # phenocam formatting
  phenocam_data <- pr_fm_phenocam(path = system.file(
    "extdata",
    package = "phenor",
    mustWork = TRUE))
  expect_type(phenocam_data, "list")

})

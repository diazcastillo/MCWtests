test_that("bMCWtest works", {
  test_temp <- tempdir()
  file.copy(system.file("extdata", "example_bMCWtest_data.csv", package = "MCWtests"), test_temp)
  bMCWtest(file.path(test_temp, "example_bMCWtest_data.csv"), 10)
  expect_true(file.exists(file.path(test_temp, "example_bMCWtest_results.csv")))
  rm(test_temp)
})

test_that("uMCWtest works", {
  test_temp <- tempdir()
  file.copy(system.file("extdata", "example_vertical_uMCWtest_data.csv", package = "MCWtests"), test_temp)
  uMCWtest(file.path(test_temp, "example_vertical_uMCWtest_data.csv"), 10)
  expect_true(file.exists(file.path(test_temp, "example_vertical_uMCWtest_results.csv")))
  rm(test_temp)
})

test_that("mbMCWtest works", {
  test_temp <- tempdir()
  file.copy(system.file("extdata", "example_vertical_mbMCWtest_data.csv", package = "MCWtests"), test_temp)
  mbMCWtest(file.path(test_temp, "example_vertical_mbMCWtest_data.csv"), 10)
  expect_true(file.exists(file.path(test_temp, "example_vertical_mbMCWtest_results.csv")))
  rm(test_temp)
})

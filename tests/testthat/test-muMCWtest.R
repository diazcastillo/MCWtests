test_that("muMCWtest works", {
  test_temp <- tempdir()
  file.copy(system.file("extdata", "example_vertical_muMCWtest_data.csv", package = "MCWtests"), test_temp)
  muMCWtest(file.path(test_temp, "example_vertical_muMCWtest_data.csv"), 10)
  expect_true(file.exists(file.path(test_temp, "example_vertical_muMCWtest_results.csv")))
  rm(test_temp)
})

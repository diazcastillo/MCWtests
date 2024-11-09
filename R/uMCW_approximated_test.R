#' Unmatched-Pairs Monte Carlo-Wilcoxon (uMCW) Approximated Test
#'
#' @param x
#' Data for uMCW approximated tests provided by function *uMCWtest*.
#'
#' @keywords internal
uMCW_approximated_test <- function(x) {

  message("running approximated tests: ")

  # running uMCW approximated tests
  uMCW_approximated_test_results <- merge(copy(x)[, uMCW_approximated_test_observed_BIs(.SD)],
                                          Reduce(function(X, Y) X[Y, on = setdiff(colnames(Y), "simulated_BI")],
                                                 pbreplicate(max_rearrangements,
                                                             copy(x)[, uMCW_approximated_test_simulated_BIs(.SD)],
                                                             simplify = FALSE)))[
                                                               , `:=`(c("expected_by_chance_BI_N", "pupper", "plower"),
                                                                      list(NCOL(.SD),
                                                                           rowSums(.SD >= observed_BI) / NCOL(.SD),
                                                                           rowSums(.SD <= observed_BI) / NCOL(.SD)))
                                                               ,.SDcols = patterns("simulated_BI")][
                                                                 , .SD, .SDcols = !patterns("simulated_BI")]
}

utils::globalVariables(codetools::findGlobals(uMCW_approximated_test, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(uMCW_approximated_test, merge = FALSE)$variables)

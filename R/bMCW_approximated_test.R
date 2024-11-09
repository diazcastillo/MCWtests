#' Biased-Measure Monte Carlo-Wilcoxon (bMCW) Approximated Test
#'
#' @param x
#' Data for bMCW approximated tests provided by function *bMCWtest*.
#'
#' @keywords internal
bMCW_approximated_test <- function(x) {

  message("running approximated tests: ")

  # running bMCW approximated tests
  bMCW_approximated_test_results <- merge(copy(x)[, bMCW_approximated_test_observed_BIs(.SD)],
                                          Reduce(function(X, Y) X[Y, on = setdiff(colnames(Y), "simulated_BI")],
                                                 pbreplicate(max_rearrangements,
                                                             copy(x)[, bMCW_approximated_test_simulated_BIs(.SD)],
                                                             simplify = FALSE)))[
                                                               , `:=`(c("expected_by_chance_BI_N", "pupper", "plower"),
                                                                      list(NCOL(.SD),
                                                                           rowSums(.SD >= observed_BI) / NCOL(.SD),
                                                                           rowSums(.SD <= observed_BI) / NCOL(.SD)))
                                                               ,.SDcols = patterns("simulated_BI")][
                                                                 , .SD, .SDcols = !patterns("simulated_BI")][
                                                                   , `:=`(subset, NULL)]
}

utils::globalVariables(codetools::findGlobals(bMCW_approximated_test, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(bMCW_approximated_test, merge = FALSE)$variables)

#' Matched-Pairs Univariate Monte Carlo-Wilcoxon (muMCW) Approximated Test
#'
#' @param x
#' Data for muMCW approximated tests provided by function *muMCWtest*.
#'
#' @keywords internal
muMCW_approximated_test <- function(x) {

  message("running approximated tests: ")

  # running muMCW approximated tests
  muMCW_approximated_test_results <- merge(copy(x)[, muMCW_approximated_test_observed_BIs(.SD)],
                                           Reduce(function(X, Y) X[Y, on = setdiff(colnames(Y), "simulated_BI")],
                                                  pbreplicate(max_rearrangements,
                                                              copy(x)[, muMCW_approximated_test_simulated_BIs(.SD)],
                                                              simplify = FALSE)))[
                                                                , `:=`(c("expected_by_chance_BI_N", "pupper", "plower"),
                                                                       list(NCOL(.SD),
                                                                            rowSums(.SD >= observed_BI) / NCOL(.SD),
                                                                            rowSums(.SD <= observed_BI) / NCOL(.SD)))
                                                                ,.SDcols = patterns("simulated_BI")][
                                                                  , .SD, .SDcols = !patterns("simulated_BI")]
}

utils::globalVariables(codetools::findGlobals(muMCW_approximated_test, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(muMCW_approximated_test, merge = FALSE)$variables)

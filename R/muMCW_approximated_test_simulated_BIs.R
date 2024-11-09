#' Matched-Pairs Univariate Monte Carlo-Wilcoxon (muMCW) Approximated Test: Simulated Bias Indexes (BIs)
#'
#' @param x
#' Data for muMCW approximated tests provided by function *muMCW_approximated_test*.
#'
#' @keywords internal
muMCW_approximated_test_simulated_BIs <- function(x) {

  # calculating BIs after randomly rearranging data
  muMCW_approximated_test_simulated_BIs <- setnames(dcast(melt(copy(x),
                                                               measure.vars = patterns("value_"),
                                                               variable.name = "pair_element",
                                                               value.name = "value")[
                                                                 , `:=`(value, value[sample(.N)]), by = internal_contrast_ID],
                                                          ... ~ pair_element,
                                                          value.var = "value")[
                                                            , muMCW_approximated_test_observed_BIs(.SD)], "observed_BI", "simulated_BI")
}

utils::globalVariables(codetools::findGlobals(muMCW_approximated_test_simulated_BIs, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(muMCW_approximated_test_simulated_BIs, merge = FALSE)$variables)

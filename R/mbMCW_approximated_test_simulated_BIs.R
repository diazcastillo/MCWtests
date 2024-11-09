#' Matched-Pairs Bivariate Monte Carlo-Wilcoxon (mbMCW) Approximated Test: Simulated Bias Indexes (BIs)
#'
#' @param x
#' Data for mbMCW approximated tests provided by function *mbMCW_approximated_test*.
#'
#' @keywords internal
mbMCW_approximated_test_simulated_BIs <- function(x) {

  # calculating BIs after randomly rearranging data
  mbMCW_approximated_test_simulated_BIs <- setnames(copy(x)[
    , `:=`(unmatched_condition, unmatched_condition[sample(.N)]), by = internal_contrast_ID][
      , mbMCW_approximated_test_observed_BIs(.SD)], "observed_BI", "simulated_BI")
}


utils::globalVariables(codetools::findGlobals(mbMCW_approximated_test_simulated_BIs, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(mbMCW_approximated_test_simulated_BIs, merge = FALSE)$variables)

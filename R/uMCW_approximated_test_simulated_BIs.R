#' Unmatched-Pairs Monte Carlo-Wilcoxon (uMCW) Approximated Test: Simulated Bias Indexes (BIs)
#'
#' @param x
#' Data for uMCW approximated tests provided by function *uMCW_approximated_test*.
#'
#' @keywords internal
uMCW_approximated_test_simulated_BIs <- function(x) {

  # calculating BIs after randomly rearranging data
  uMCW_approximated_test_simulated_BIs <- setnames(copy(x)[
    , `:=`(value, value[sample(.N)]), by = internal_contrast_ID][
      , uMCW_approximated_test_observed_BIs(.SD)], "observed_BI", "simulated_BI")
}

utils::globalVariables(codetools::findGlobals(uMCW_approximated_test_simulated_BIs, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(uMCW_approximated_test_simulated_BIs, merge = FALSE)$variables)

#' Biased-Measure Monte Carlo-Wilcoxon (bMCW) Approximated Test: Simulated Bias Indexes (BIs)
#'
#' @param x
#' Data for bMCW approximated tests provided by function *bMCW_approximated_test*.
#'
#' @keywords internal
bMCW_approximated_test_simulated_BIs <- function(x) {

  # calculating BIs after randomly rearranging data
  bMCW_approximated_test_simulated_BIs <- setnames(rbind(if (NROW(unique(copy(x)[subset_type == "wholeset", internal_contrast_ID])) > 0) {
    copy(x)[subset_type == "wholeset", ][
      , `:=`(sign_bias_value, sign_bias_value[sample(.N)]), by = internal_contrast_ID][
        , bMCW_approximated_test_observed_BIs(.SD)]
  }, if (NROW(unique(copy(x)[subset_type != "wholeset", internal_contrast_ID])) > 0) {
    copy(x)[subset_type != "wholeset", ][
      , `:=`(subset, subset[sample(.N)]), by = internal_contrast_ID][
        , bMCW_approximated_test_observed_BIs(.SD)]
  }), "observed_BI", "simulated_BI")
}

utils::globalVariables(codetools::findGlobals(bMCW_approximated_test_simulated_BIs, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(bMCW_approximated_test_simulated_BIs, merge = FALSE)$variables)

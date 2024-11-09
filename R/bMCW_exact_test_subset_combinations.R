#' Biased-Measure Monte Carlo-Wilcoxon (bMCW) Exact Test: All Subset Combinations
#'
#' @param x
#' Data for bMCW exact tests provided by function *bMCW_exact_test*.
#'
#' @keywords internal
bMCW_exact_test_subset_combinations <- function(x) {

  # drawing all possible data rearrangements for subset bMCW exact tests
  bMCW_exact_test_subset_combinations <- setorder(melt(data.table(t(combn(copy(x)[, internal_element_ID], NROW(copy(x)[subset != "NO", ]))))[
    , `:=`(combination = paste0("C", seq_len(NROW(V1))))],
    measure.vars = patterns("V"),
    variable.name = "x",
    value.name = "internal_element_ID")[
      , `:=`(x, NULL)],
    combination)
}

utils::globalVariables(codetools::findGlobals(bMCW_exact_test_subset_combinations, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(bMCW_exact_test_subset_combinations, merge = FALSE)$variables)

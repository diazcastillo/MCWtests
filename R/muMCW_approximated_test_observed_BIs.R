#' Matched-Pairs Univariate Monte Carlo-Wilcoxon (muMCW) Approximated Test: Observed Bias Indexes (BIs)
#'
#' @param x
#' Data for muMCW approximated tests provided by function *muMCW_approximated_test*.
#'
#' @keywords internal
muMCW_approximated_test_observed_BIs <- function(x) {

  # formatting muMCW approximated test dataset
  muMCW_approximated_test_signed_ranks <- melt(copy(x)[
    , `:=`(`a-b_condition_contrast` = paste0(condition_a, "-", condition_b),
           `b-a_condition_contrast` = paste0(condition_b, "-", condition_a),
           abs_diff = abs(value_a - value_b))][
             abs_diff == 0, `:=`(abs_diff, NA)][
               , `:=`(abs_rank = rank(abs_diff, na.last = "keep", ties.method = "min")),
               by = .(internal_contrast_ID)][
                 , `:=`(`a-b_signed_rank` = abs_rank * sign(value_a - value_b),
                        `b-a_signed_rank` = abs_rank * sign(value_b - value_a))][
                          is.na(`a-b_signed_rank`), `:=`(`a-b_signed_rank`, 0)][
                            is.na(`b-a_signed_rank`), `:=`(`b-a_signed_rank`, 0)],
    id.vars = "internal_contrast_ID",
    measure.vars = patterns("condition_contrast", "_signed_rank$"),
    variable.name = "condition",
    value.name = c("condition_contrast", "signed_rank"))

  # calculating observed BIs
  muMCW_approximated_test_observed_BIs <- copy(muMCW_approximated_test_signed_ranks)[
    , .(BI_type = "muMCW_BI", N = .N, sum_signed_ranks = sum(signed_rank)), by = .(internal_contrast_ID, condition_contrast)][
      , `:=`(observed_BI = (sum_signed_ranks / ((N * (N + 1)) / 2)))][
        , `:=`(sum_signed_ranks = NULL)][
          , `:=`(N, NULL)]

  return (muMCW_approximated_test_observed_BIs)
}

utils::globalVariables(codetools::findGlobals(muMCW_approximated_test_observed_BIs, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(muMCW_approximated_test_observed_BIs, merge = FALSE)$variables)

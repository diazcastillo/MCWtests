#' Matched-Pairs Bivariate Monte Carlo-Wilcoxon (mbMCW) Approximated Test: Observed Bias Indexes (BIs)
#'
#' @param x
#' Data for mbMCW approximated tests provided by function *mbMCW_approximated_test*.
#'
#' @keywords internal
mbMCW_approximated_test_observed_BIs <- function(x) {

  # calculating observed BIs
  mbMCW_approximated_test_observed_BIs <- melt(dcast(melt(copy(x)[, .(sum_abs_rank = sum(abs_rank)), by = internal_contrast_ID][
    copy(x)[, .(`a-b_sum_signed_rank` = sum(`a-b_signed_rank`), `b-a_sum_signed_rank` = sum(`b-a_signed_rank`)), by = .(internal_contrast_ID, unmatched_condition)],
    on = .(internal_contrast_ID)],
    id.vars = c("internal_contrast_ID", "unmatched_condition", "sum_abs_rank"),
    measure.vars = patterns("_signed_rank"),
    variable.name = "matched_condition_contrast",
    value.name = "matched_sum_signed_ranks")[
      , `:=`(matched_condition_contrast, gsub("_sum_signed_rank", "", matched_condition_contrast))],
    internal_contrast_ID + matched_condition_contrast + sum_abs_rank ~ unmatched_condition,
    value.var = c("matched_sum_signed_ranks"))[
      , `:=`(`x-y_BI` = (x - y) / sum_abs_rank,
             `y-x_BI` = (y - x) / sum_abs_rank)],
    id.vars = c("internal_contrast_ID", "matched_condition_contrast"),
    measure.vars = patterns("_BI"),
    variable.name = "unmatched_condition_contrast",
    value.name = "observed_BI")[
      , `:=`(unmatched_condition_contrast, gsub("_BI", "", unmatched_condition_contrast))]
  }

utils::globalVariables(codetools::findGlobals(mbMCW_approximated_test_observed_BIs, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(mbMCW_approximated_test_observed_BIs, merge = FALSE)$variables)

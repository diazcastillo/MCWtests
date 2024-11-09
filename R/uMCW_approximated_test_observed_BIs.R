#' Unmatched-Pairs Monte Carlo-Wilcoxon (uMCW) Approximated Test: Observed Bias Indexes (BIs)
#'
#' @param x
#' Data for uMCW approximated tests provided by function *uMCW_approximated_test*.
#'
#' @keywords internal
uMCW_approximated_test_observed_BIs <- function(x) {

  # formatting uMCW approximated test dataset
  uMCW_approximated_test_signed_ranks <- melt(copy(x)[, `:=`(aid, paste0(replicate, "/", condition, "/", value))][
    , CJ(aid_a = aid, aid_b = aid), by = internal_contrast_ID][
      , `:=`(c("replicate_a", "condition_a", "value_a"), tstrsplit(aid_a, "/", fixed = TRUE))][
        , `:=`(c("replicate_b", "condition_b", "value_b"), tstrsplit(aid_b, "/", fixed = TRUE))][
          , `:=`(value_a = as.numeric(value_a), value_b = as.numeric(value_b))][
            replicate_a != replicate_b, ][
              , `:=`(aid_a = NULL, aid_b = NULL)][
                , `:=`(BI = "uMCW_BI", HBI = "uMCW_HBI")],
    measure.vars = patterns("BI"),
    value.name = "BI_type")[
      , `:=`(variable = NULL)][
        BI_type == "uMCW_BI" | (BI_type == "uMCW_HBI" & condition_a == condition_b), ][
          , `:=`(aid = fifelse(BI_type == "uMCW_BI",
                               paste0(internal_contrast_ID, ".", BI_type, ".", replicate_a, ".", replicate_b),
                               paste0(internal_contrast_ID, ".", BI_type, ".", pmin(replicate_a, replicate_b), ".", pmax(replicate_a, replicate_b))))][
                                 !duplicated(aid)][
                                   , `:=`(abs_diff = abs(value_a - value_b))][
                                     abs_diff == 0, `:=`(abs_diff = NA)][
                                       , `:=`(rank = rank(abs_diff, na.last = "keep", ties.method = "min")),
                                       by = .(internal_contrast_ID, BI_type)][
                                         is.na(rank), `:=`(rank = 0)][
                                           , `:=`(signed_rank = fifelse(BI_type == "uMCW_BI", rank * sign(value_a - value_b), rank))]

  # calculating observed BIs
  uMCW_approximated_test_observed_BIs <- setorder(rbind(merge(copy(uMCW_approximated_test_signed_ranks)[BI_type == "uMCW_BI", .(N = .N), by = .(internal_contrast_ID)],
                                                              copy(uMCW_approximated_test_signed_ranks)[BI_type == "uMCW_BI" & condition_a != condition_b, .(condition_contrast = paste0(condition_a, "-", condition_b),
                                                                                                                                                                n = .N,
                                                                                                                                                                sum_rank = sum(signed_rank)),
                                                                                                        by = .(internal_contrast_ID, BI_type, condition_a, condition_b)])[
                                                                                                          , `:=`(observed_BI = (sum_rank / (((N * (N + 1)) / 2) - (((N - n) * ((N - n) + 1)) / 2))), condition_a = NULL, condition_b = NULL, sum_rank = NULL, N = NULL, n = NULL)],
                                                        melt(dcast(merge(copy(uMCW_approximated_test_signed_ranks)[BI_type == "uMCW_HBI", .(sum_abs_rank = sum(rank)), by = .(internal_contrast_ID)],
                                                                         copy(uMCW_approximated_test_signed_ranks)[BI_type == "uMCW_HBI", .(sum_abs_rank_per_condition = sum(rank)),
                                                                                                                   by = .(internal_contrast_ID, BI_type, condition_a, condition_b)])[
                                                                                                                     , `:=`(condition = condition_a, condition_a = NULL, condition_b = NULL)],
                                                                   internal_contrast_ID + BI_type + sum_abs_rank ~ condition,
                                                                   value.var = "sum_abs_rank_per_condition")[
                                                                     , `:=`(`a-b_condition` = "a-b",
                                                                            `a-b_HBI` = (a - b) / sum_abs_rank,
                                                                            `b-a_condition` = "b-a",
                                                                            `b-a_HBI` = (b - a) / sum_abs_rank)],
                                                             id.vars = c("internal_contrast_ID", "BI_type"),
                                                             measure.vars = patterns("_condition", "_HBI"),
                                                             value.name = c("condition_contrast", "observed_BI"))[
                                                               , `:=`(variable = NULL)]),
                                                  internal_contrast_ID, BI_type, condition_contrast)

  return (uMCW_approximated_test_observed_BIs)
}

utils::globalVariables(codetools::findGlobals(uMCW_approximated_test_observed_BIs, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(uMCW_approximated_test_observed_BIs, merge = FALSE)$variables)

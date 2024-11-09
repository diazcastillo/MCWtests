#' Biased-Measure Monte Carlo-Wilcoxon (bMCW) Exact Test
#'
#' @param x
#' Data for bMCW exact tests provided by function *bMCWtest*.
#'
#' @keywords internal
bMCW_exact_test <- function(x) {

  message("running exact tests: ")

  # drawing all possible data rearrangements
  bMCW_exact_test_combinations <- setorder(rbind(
    if (NROW(unique(copy(x)[subset_type == "wholeset", internal_contrast_ID])) > 0) {
      rbind(copy(x)[subset_type == "wholeset", .(internal_contrast_ID, sign_bias_value, rank_abs_bias_value)][
        , `:=`(combination = "C0")],
        copy(x)[subset_type == "wholeset", bMCW_exact_test_wholeset_combinations(.SD), by = internal_contrast_ID][
          , `:=`(replicate_pair = NULL)])[
            , `:=`(subset_type = "wholeset", signed_rank = rank_abs_bias_value * sign_bias_value, rank_abs_bias_value = NULL, sign_bias_value = NULL)]
    }, if (NROW(unique(copy(x)[subset_type != "wholeset", internal_contrast_ID])) > 0) {
      rbind(copy(x)[subset_type != "wholeset" & subset != "NO", .(internal_contrast_ID, internal_element_ID, combination = "C0")],
            copy(x)[subset_type != "wholeset", ][
              , bMCW_exact_test_subset_combinations(.SD), by = internal_contrast_ID])[
                copy(x)[subset_type != "wholeset"][
                  , `:=`(signed_rank = rank_abs_bias_value * sign_bias_value)][
                    , .(internal_contrast_ID, internal_element_ID, signed_rank)]
                , on = .(internal_contrast_ID = internal_contrast_ID, internal_element_ID = internal_element_ID)][
                  , `:=`(subset_type = "subset", internal_element_ID = NULL)]
    }), internal_contrast_ID, combination)

  # calculating BIs
  bMCW_exact_test_BI <- setorder(rbind(copy(bMCW_exact_test_combinations)[subset_type == "wholeset", .(N = .N, sum_signed_ranks = sum(signed_rank)), by = .(internal_contrast_ID, combination)][
    , `:=`(observed_BI = (sum_signed_ranks / ((N * (N + 1)) / 2)))][
      , `:=`(BI_type = "bMCW_wBI", N = NULL, sum_signed_ranks = NULL)],
    copy(bMCW_exact_test_combinations)[subset_type != "wholeset", .(n = .N, sum_signed_ranks = sum(signed_rank)), by = .(internal_contrast_ID, combination)][
      unique(copy(x)[subset_type != "wholeset", ][
        , .(internal_contrast_ID, N)]), on = .(internal_contrast_ID = internal_contrast_ID)][
          , `:=`(observed_BI = sum_signed_ranks / (((N * (N + 1)) / 2) - (((N - n) * ((N - n) + 1)) / 2)))][
            , `:=`(BI_type = "bMCW_sBI", N = NULL, n = NULL, sum_signed_ranks = NULL)]), internal_contrast_ID, combination)

  # writing bMCW exact test result table
  bMCW_exact_test_results <- setorder(rbindlist(lapply(split(copy(bMCW_exact_test_BI), bMCW_exact_test_BI$internal_contrast_ID),
                                                       function(x) {
                                                         a <- dcast(x, internal_contrast_ID + BI_type ~ combination, value.var = "observed_BI")[
                                                           , `:=`(c("observed_BI", "expected_by_chance_BI_N", "pupper", "plower"),
                                                                  list(C0,
                                                                       NCOL(.SD),
                                                                       rowSums(.SD >= C0) / NCOL(.SD),
                                                                       rowSums(.SD <= C0) / NCOL(.SD))),
                                                           .SDcols = patterns("^C([1-9])")][
                                                             , .SD, .SDcols = !patterns ("^C")]
                                                       })), internal_contrast_ID)

  return(bMCW_exact_test_results)
}

utils::globalVariables(codetools::findGlobals(bMCW_exact_test, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(bMCW_exact_test, merge = FALSE)$variables)

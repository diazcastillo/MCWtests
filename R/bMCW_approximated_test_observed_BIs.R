#' Biased-Measure Monte Carlo-Wilcoxon (bMCW) Approximated Test: Observed Bias Indexes (BIs)
#'
#' @param x
#' Data for bMCW approximated tests provided by function *bMCW_approximated_test*.
#'
#' @keywords internal
bMCW_approximated_test_observed_BIs <- function(x) {

  # calculating observed BIs
  bMCW_approximated_test_observed_BIs <- rbind(
    if(NROW(unique(copy(x)[subset_type == "wholeset", internal_contrast_ID])) > 0) {
      copy(x)[subset_type == "wholeset", ][
        ,`:=`(signed_rank = rank_abs_bias_value * sign_bias_value)][
          , .(N = .N, sum_signed_ranks = sum(signed_rank)), by = internal_contrast_ID][
            , `:=`(observed_BI = (sum_signed_ranks / ((N * (N + 1)) / 2)))][
              , `:=`(subset = "none", BI_type = "bMCW_wBI", N = NULL, sum_signed_ranks = NULL)]
    }, if(NROW(unique(copy(x)[subset_type != "wholeset", internal_contrast_ID])) > 0) {
      merge(copy(x)[subset_type != "wholeset", ][
        ,`:=`(signed_rank = rank_abs_bias_value * sign_bias_value)][
          , .(N = .N), by = internal_contrast_ID],
        copy(x)[subset_type != "wholeset", ][
          ,`:=`(signed_rank = rank_abs_bias_value * sign_bias_value)][
            , .(n = .N, sum_signed_ranks = sum(signed_rank)), by = .(internal_contrast_ID, subset)])[
              subset != "NO", ][
                , `:=`(observed_BI, (sum_signed_ranks / (((N * (N + 1)) / 2) - (((N - n) * ((N - n) + 1)) / 2))))][
                  , `:=`(BI_type = "bMCW_sBI", N = NULL, n = NULL, sum_signed_ranks = NULL)]
    })
}

utils::globalVariables(codetools::findGlobals(bMCW_approximated_test_observed_BIs, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(bMCW_approximated_test_observed_BIs, merge = FALSE)$variables)

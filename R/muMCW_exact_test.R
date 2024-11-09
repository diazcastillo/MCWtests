#' Matched-Pairs Univariate Monte Carlo-Wilcoxon (muMCW) Exact Test
#'
#' @param x
#' Data for muMCW exact tests provided by function *muMCWtest*.
#'
#' @keywords internal
muMCW_exact_test <- function(x) {

  message("running exact tests: ")

  # setting muMCW_exact_test dataset
  muMCW_exact_test_dataset <- copy(x)[
    , `:=`(pair_element_a = as.character(seq(from = 1, to = .N * 2, by = 2)),
           pair_element_b = as.character(seq(from = 2, to = .N * 2, by = 2)))
    , by = .(internal_contrast_ID)]

  # drawing all possible data rearrangements
  muMCW_exact_test_combinations <- melt(rbind(copy(muMCW_exact_test_dataset)[
    , .(internal_contrast_ID, condition_a, condition_b, replicate_pair, value_a, value_b)][
      , `:=`(combination = "C0")],
    dcast(merge(unique(copy(muMCW_exact_test_dataset)[, .(internal_contrast_ID, N)])[
      melt(setnames(data.table(N = unique(copy(muMCW_exact_test_dataset)[, N]))[
        , `:=`(aid, N)][
          , muMCW_exact_test_combinations(.SD), by = aid], "aid", "N"),
        measure.vars = patterns("pair_element_"),
        variable.name = "pair_element",
        value.name = "replicate_ID"), on = .(N = N), allow.cartesian = TRUE],
      melt(copy(muMCW_exact_test_dataset),
           measure.var = patterns("pair_element_", "value_"),
           value.name = c("replicate_ID", "value"))[
             , .(internal_contrast_ID, condition_a, condition_b, replicate_ID, value)]),
      internal_contrast_ID + combination + condition_a + condition_b + replicate_pair ~ pair_element, value.var = "value")[
        , `:=`(value_a = pair_element_a, value_b = pair_element_b, pair_element_a = NULL, pair_element_b = NULL)])[
          , `:=`(`a-b_condition_contrast` = paste0(condition_a, "-", condition_b),
                 `b-a_condition_contrast` = paste0(condition_b, "-", condition_a),
                 abs_diff = abs(value_a - value_b))][
                   abs_diff == 0, `:=`(abs_diff, NA)][
                     , `:=`(abs_rank = rank(abs_diff, na.last = "keep", ties.method = "min")),
                     by = .(internal_contrast_ID, combination)][
                       , `:=`(`a-b_signed_rank` = abs_rank * sign(value_a - value_b),
                              `b-a_signed_rank` = abs_rank * sign(value_b - value_a))][
                                is.na(`a-b_signed_rank`), `:=`(`a-b_signed_rank`, 0)][
                                  is.na(`b-a_signed_rank`), `:=`(`b-a_signed_rank`, 0)],
    id.vars = c("internal_contrast_ID", "combination"),
    measure.vars = patterns("condition_contrast", "_signed_rank$"),
    variable.name = "condition",
    value.name = c("condition_contrast", "signed_rank"))

  # calculating BIs
  muMCW_exact_test_BI <- copy(muMCW_exact_test_combinations)[
    , .(BI_type = "muMCW_BI", N = .N, sum_signed_ranks = sum(signed_rank)), by = .(internal_contrast_ID, combination, condition_contrast)][
      , `:=`(observed_BI = (sum_signed_ranks / ((N * (N +1)) / 2)))]

  # writing muMCW exact test result table
  muMCW_exact_test_results <- setorder(rbindlist(lapply(split(copy(muMCW_exact_test_BI), muMCW_exact_test_BI$internal_contrast_ID),
                                                        function(x) {
                                                          a <- dcast(x, internal_contrast_ID + BI_type + condition_contrast ~ combination, value.var = "observed_BI")[
                                                            , `:=`(c("observed_BI", "expected_by_chance_BI_N", "pupper", "plower"),
                                                                   list(C0,
                                                                        NCOL(.SD),
                                                                        rowSums(.SD >= C0) / NCOL(.SD),
                                                                        rowSums(.SD <= C0) / NCOL(.SD))),
                                                            .SDcols = patterns("^C([1-9])")][
                                                              , .SD, .SDcols = !patterns ("^C")]
                                                        })), internal_contrast_ID, condition_contrast)

  return(muMCW_exact_test_results)
}

utils::globalVariables(codetools::findGlobals(muMCW_exact_test, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(muMCW_exact_test, merge = FALSE)$variables)

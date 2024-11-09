#' Matched-Pairs Bivariate Monte Carlo-Wilcoxon (mbMCW) Exact Test
#'
#' @param x
#' Data for mbMCW exact tests provided by function *mbMCWtest*.
#'
#' @keywords internal
mbMCW_exact_test <- function(x) {

  message("running exact tests: ")

  # drawing all possible data rearrangements
  mbMCW_exact_test_combinations <- rbind(copy(x)[, .(combination = "C0", internal_contrast_ID, replicate_pair, abs_rank, unmatched_condition, `a-b_condition_contrast`, `b-a_condition_contrast`, `a-b_signed_rank`, `b-a_signed_rank`)][
    , `:=`(replicate_pair, as.numeric(replicate_pair))][
      , `:=`(unmatched_condition, fifelse(unmatched_condition == unique(unmatched_condition)[1], "x", "y")), by = internal_contrast_ID],
    merge(copy(x)[, CJ(combination = paste0("C", seq_len(ncol(combn(unique(N), unique(N_x))))),
                       replicate_pair = paste0(replicate_pair, "/", unmatched_condition)), by = internal_contrast_ID][
                         , `:=`(c("replicate_pair", "reference_condition"), tstrsplit(replicate_pair, "/", fixed = TRUE))][
                           , `:=`(replicate_pair, as.numeric(replicate_pair))],
          na.omit(melt(copy(x)[, combn(replicate_pair, unique(N_x), FUN = paste, collapse = "/"), by = internal_contrast_ID][
            , `:=`(combination = paste0("C", seq_len(.N)), unmatched_condition = "x"), by = internal_contrast_ID][
              , `:=`(paste0("V", seq_len(max(copy(x)[, N_x]))), tstrsplit(V1, "/"))],
            measure.vars = patterns("V"),
            value.name = "replicate_pair"))[
              , `:=`(replicate_pair = as.numeric(replicate_pair), variable = NULL)],
          all.x = TRUE, allow.cartesian = TRUE)[
            is.na(unmatched_condition), `:=`(unmatched_condition = "y")][
              , `:=`(reference_condition, NULL)][
                copy(x)[, .(internal_contrast_ID, replicate_pair, abs_rank, `a-b_condition_contrast`, `b-a_condition_contrast`, `a-b_signed_rank`, `b-a_signed_rank`)][
                  , `:=`(replicate_pair, as.numeric(replicate_pair))], on = .(internal_contrast_ID, replicate_pair), allow.cartesian = TRUE])

  # calculating BIs
  mbMCW_exact_test_BI <- melt(dcast(melt(copy(mbMCW_exact_test_combinations)[, .(sum_abs_rank = sum(abs_rank)), by = .(internal_contrast_ID, combination)][
    copy(mbMCW_exact_test_combinations)[, .(`a-b_sum_signed_rank` = sum(`a-b_signed_rank`), `b-a_sum_signed_rank` = sum(`b-a_signed_rank`)), by = .(internal_contrast_ID, combination, unmatched_condition)],
    on = .(internal_contrast_ID, combination)],
    id.vars = c("internal_contrast_ID", "combination", "unmatched_condition", "sum_abs_rank"),
    measure.vars = patterns("_signed_rank"),
    variable.name = "matched_condition_contrast",
    value.name = "matched_sum_signed_ranks")[
      , `:=`(matched_condition_contrast, gsub("_sum_signed_rank", "", matched_condition_contrast))],
    internal_contrast_ID + combination + matched_condition_contrast + sum_abs_rank ~ unmatched_condition,
    value.var = c("matched_sum_signed_ranks"))[
      , `:=`(`x-y_BI` = (x - y) / sum_abs_rank,
             `y-x_BI` = (y - x) / sum_abs_rank)],
    id.vars = c("internal_contrast_ID", "combination", "matched_condition_contrast"),
    measure.vars = patterns("_BI"),
    variable.name = "unmatched_condition_contrast",
    value.name = "observed_BI")[
      , `:=`(unmatched_condition_contrast, gsub("_BI", "", unmatched_condition_contrast))][
        , `:=`(BI_type = "mbMCW_BI")]

  # writing mbMCW exact test result table
  mbMCW_exact_test_results <- setorder(rbindlist(lapply(split(copy(mbMCW_exact_test_BI), mbMCW_exact_test_BI$internal_contrast_ID),
                                                        function(x) {
                                                          a <- dcast(x, internal_contrast_ID + matched_condition_contrast + unmatched_condition_contrast + BI_type ~ combination, value.var = "observed_BI")[
                                                            , `:=`(c("expected_by_chance_BI_N", "observed_BI", "pupper", "plower"),
                                                                   list(NCOL(.SD),
                                                                        C0,
                                                                        rowSums(.SD >= C0) / NCOL(.SD),
                                                                        rowSums(.SD <= C0) / NCOL(.SD))),
                                                            .SDcols = patterns("^C([1-9])")][
                                                              , .SD, .SDcols = !patterns ("^C")]
                                                        })), internal_contrast_ID, matched_condition_contrast, unmatched_condition_contrast)

  return(mbMCW_exact_test_results)
}

utils::globalVariables(codetools::findGlobals(mbMCW_exact_test, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(mbMCW_exact_test, merge = FALSE)$variables)

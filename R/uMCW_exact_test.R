#' Unmatched-Pairs Monte Carlo-Wilcoxon (uMCW) Exact Test
#'
#' @param x
#' Data for uMCW exact tests provided by function *uMCWtest*.
#'
#' @keywords internal
uMCW_exact_test <- function(x) {

  message("running exact tests: ")

  # drawing all possible data rearrangements
  uMCW_exact_test_combinations <- melt(rbind(copy(x)[, `:=`(aid, paste0(replicate, "/", condition, "/", value))][
    , CJ(aid_a = aid, aid_b = aid), by = internal_contrast_ID][
      , `:=`(combination = "C0")],
    merge(copy(x)[, CJ(combination = paste0("C", seq_len(ncol(combn(unique(N), unique(n_a))))),
                       replicate = paste0(replicate, "/", condition)), by = internal_contrast_ID][
                         , `:=`(c("replicate", "reference_condition"), tstrsplit(replicate, "/", fixed = TRUE))][
                           , `:=`(replicate, as.numeric(replicate))],
          na.omit(melt(copy(x)[, combn(replicate, unique(n_a), FUN = paste, collapse = "/"), by = internal_contrast_ID][
            , `:=`(combination = paste0("C", seq_len(.N)), condition = "a"), by = internal_contrast_ID][
              , paste0("V", seq_len(max(copy(x)[, n_a]))) := tstrsplit(V1, "/")],
            measure.vars = patterns("V"),
            value.name = "replicate"))[
              , `:=`(replicate = as.numeric(replicate), variable = NULL)],
          all.x = TRUE, allow.cartesian = TRUE)[
            , `:=`(condition, fcase(condition == "a", unique(reference_condition)[1],
                                    is.na(condition), unique(reference_condition)[2])), by = internal_contrast_ID][
                                      , `:=`(reference_condition, NULL)][
                                        copy(x)[, .(internal_contrast_ID, replicate, value)], on = .(internal_contrast_ID, replicate), allow.cartesian = TRUE][
                                          , `:=`(aid, paste0(replicate, "/", condition, "/", value))][
                                            , CJ(aid_a = aid, aid_b = aid), by = .(internal_contrast_ID, combination)])[
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
                               paste0(internal_contrast_ID, combination, ".", BI_type, ".", replicate_a, ".", replicate_b),
                               paste0(internal_contrast_ID, combination, ".", BI_type, ".", pmin(replicate_a, replicate_b), ".", pmax(replicate_a, replicate_b))))][
                                 !duplicated(aid)][
                                   , `:=`(abs_diff = abs(value_a - value_b))][
                                     abs_diff == 0, `:=`(abs_diff = NA)][
                                       , `:=`(rank = rank(abs_diff, na.last = "keep", ties.method = "min")),
                                       by = .(internal_contrast_ID, combination, BI_type)][
                                         is.na(rank), `:=`(rank = 0)][
                                           , `:=`(signed_rank = fifelse(BI_type == "uMCW_BI", rank * sign(value_a - value_b), rank))]

  # calculating BIs
  uMCW_exact_test_BI <- setorder(rbind(merge(copy(uMCW_exact_test_combinations)[BI_type == "uMCW_BI", .(N = .N), by = .(internal_contrast_ID, combination)],
                                             copy(uMCW_exact_test_combinations)[BI_type == "uMCW_BI" & condition_a != condition_b, .(condition_contrast = paste0(condition_a, "-", condition_b),
                                                                                                                                        n = .N,
                                                                                                                                        sum_rank = sum(signed_rank)),
                                                                                by = .(internal_contrast_ID, combination, BI_type, condition_a, condition_b)])[
                                                                                  , `:=`(observed_BI = (sum_rank / (((N * (N + 1)) / 2) - (((N - n) * ((N - n) + 1)) / 2))), condition_a = NULL, condition_b = NULL, sum_rank = NULL, N = NULL, n = NULL)],
                                       melt(dcast(merge(copy(uMCW_exact_test_combinations)[BI_type == "uMCW_HBI", .(sum_abs_rank = sum(rank)), by = .(internal_contrast_ID, combination)],
                                                        copy(uMCW_exact_test_combinations)[BI_type == "uMCW_HBI", .(sum_abs_rank_per_condition = sum(rank)),
                                                                                           by = .(internal_contrast_ID, combination, BI_type, condition_a, condition_b)])[
                                                                                             , `:=`(condition = condition_a, condition_a = NULL, condition_b = NULL)],
                                                  internal_contrast_ID + combination + BI_type + sum_abs_rank ~ condition,
                                                  value.var = "sum_abs_rank_per_condition")[
                                                    , `:=`(`a-b_condition` = "a-b",
                                                           `a-b_HBI` = (a - b) / sum_abs_rank,
                                                           `b-a_condition` = "b-a",
                                                           `b-a_HBI` = (b - a) / sum_abs_rank)],
                                            id.vars = c("internal_contrast_ID", "combination", "BI_type"),
                                            measure.vars = patterns("_condition", "_HBI"),
                                            value.name = c("condition_contrast", "observed_BI"))[
                                              , `:=`(variable = NULL)]), internal_contrast_ID, combination, BI_type, condition_contrast)

  # writing uMCW exact test result table
  uMCW_exact_test_results <- setorder(rbindlist(lapply(split(copy(uMCW_exact_test_BI), uMCW_exact_test_BI$internal_contrast_ID),
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

  return (uMCW_exact_test_results)
}

utils::globalVariables(codetools::findGlobals(uMCW_exact_test, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(uMCW_exact_test, merge = FALSE)$variables)

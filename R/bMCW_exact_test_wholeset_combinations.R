#' Biased-Measure Monte Carlo-Wilcoxon (bMCW) Exact Test: All Whole-set Combinations
#'
#' @param x
#' Data for bMCW exact tests provided by function *bMCW_exact_test*.
#'
#' @keywords internal
bMCW_exact_test_wholeset_combinations <- function(x) {

  # drawing all possible data rearrangements for wholeset bMCW exact tests
  bMCW_exact_test_wholeset_combinations <- dcast(melt(melt(data.table(seed = copy(x)[
    , CJ(unique(copy(x)[, element_rank]),
         unique(copy(x)[, element_sign]))][
           , `:=`(pair = paste0(V1, ":", V2), V1 = NULL, V2 = NULL)][
             , combn(pair, NROW(copy(x)), paste, collapse = "/")])[
               , `:=`(V1, gsub("/", ":", seed))][
                 , tstrsplit(V1, ":", fixed = TRUE), by = seed][
                   , `:=`(count, apply(.SD, 1, uniqueN)), .SDcols = patterns("V")][
                     count == unique(copy(x)[, N]) * 2, ][
                       , tstrsplit(seed,"/", fixed = TRUE)][
                         , `:=`(combination, paste0("C", .I))],
    measure.vars = patterns("^V"),
    variable.name = "replicate_pair",
    value.name = "pair")[
      , `:=`(c("element_rank", "element_sign"), tstrsplit(pair, ":", fixed = TRUE))],
    id.vars = c("combination", "replicate_pair"),
    measure.vars = patterns("element"),
    variable.name = "element_trait_name",
    value.name = "element_trait")[
      melt(copy(x)[, .(element_rank, element_sign, sign_bias_value, rank_abs_bias_value)],
           measure.vars = patterns("element", "bias_value"),
           value.name = c("element_trait", "element_trait_value")), on = .(element_trait = element_trait)],
    combination + replicate_pair ~ element_trait_name,
    value.var = "element_trait_value")[
      , `:=`(sign_bias_value = element_sign, rank_abs_bias_value = element_rank, element_sign = NULL, element_rank = NULL)]
}

utils::globalVariables(codetools::findGlobals(bMCW_exact_test_wholeset_combinations, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(bMCW_exact_test_wholeset_combinations, merge = FALSE)$variables)

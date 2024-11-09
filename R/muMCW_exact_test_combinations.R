#' Matched-Pairs Univariate Monte Carlo-Wilcoxon (muMCW) Exact Test: Data Combinations
#'
#' @param x
#' Data for muMCW exact tests provided by function *muMCW_exact_test*.
#'
#' @keywords internal
muMCW_exact_test_combinations <- function(x) {

  # drawing all possible data rearrangements
  muMCW_exact_test_combinations <- setorder(melt(data.table(seed = combn(data.table(t(combn(1:(copy(x)[, N] * 2), 2, paste, collapse = ":")))[1], copy(x)[, N], paste, collapse = "/"))[
    , `:=`(V1, gsub("/", ":", seed))][
      , tstrsplit(V1, ":", fixed = TRUE), by = seed][
        , `:=`(count, apply(.SD, 1, uniqueN)), .SDcols = patterns("V")][
          count == copy(x)[, N] * 2, ][
            , tstrsplit(seed,"/", fixed=TRUE)][
              ,`:=`(combination = paste0("C", .I))],
    id.vars = "combination", variable.name = "replicate_pair", value.name = "pair")[
      , `:=`(c("pair_element_a", "pair_element_b"), tstrsplit(pair, ":", fixed = TRUE))][
        , `:=`(pair, NULL)][
          , `:=`(replicate_pair, gsub("V", "", replicate_pair))], combination, replicate_pair)
}

utils::globalVariables(codetools::findGlobals(muMCW_exact_test_combinations, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(muMCW_exact_test_combinations, merge = FALSE)$variables)

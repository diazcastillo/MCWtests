#' @title
#' Matched-Measures Univariate Monte Carlo-Wilcoxon (muMCW) Test
#'
#' @description
#' The muMCW test is a statistical tool to assess whether one set of inherently matched-paired measures is significantly biased in the same direction. For instance, muMCW tests can be used to analyze bodyweights or transcript abundances determined at two different timepoints for the same set of mice.
#'
#' @details
#' The function *muMCWtest* eliminates any matched-paired measures with at least one missing value (NA) before proceeding with the following steps.
#'
#' - To estimate the bias for all matched-paired measures in the dataset, the function *muMCWtest* performs the following tasks:
#'
#'    - For each matched-pair of measures, it subtracts values for the two possible condition contrasts (*e.g.*, *a-b* and *b-a*).
#'    - For each condition contrast, it ranks the absolute values of non-zero differences from lowest to highest. Measure pair differences with a value of 0 are assigned a 0 rank. If multiple measure pair differences have the same absolute value, all tied measure pair differences are assigned the lowest rank possible.
#'    - It assigns each measure pair rank a sign based on the sign of its corresponding measure pair difference.
#'    - It sums the signed ranks for each condition contrast.
#'    - It calculates muMCW_BI by dividing each sum of signed ranks by the maximum number that sum could have if the corresponding measure pairs had the highest possible positive ranks. Consequently, muMCW_BI ranges between 1 when all measures corresponding to the first condition are higher than all measures corresponding to the second condition, and -1 when all measures corresponding to the first condition are lower than all measures corresponding to the second condition.
#'
#'  - To assess the significance of the muMCW_BIs obtained from the user-provided dataset (observed muMCW_BIs), the function *muMCWtest* performs the following tasks:
#'
#'    - It generates a collection of expected-by-chance muMCW_BIs. These expected values are obtained by rearranging the measures between and within the two conditions multiple times. The user-provided parameter *max_rearrangements* determines the two paths that the function *muMCWtest* can follow to generate the collection of expected-by-chance muMCW_BIs:
#'
#'      - *muMCW exact testing*: If the number of distinct measure rearrangements that can alter their initial pair and set distribution is less than *max_rearrangements*, the function *muMCWtest* calculates muMCW_BIs for all possible data rearrangements.
#'      - *muMCW approximated testing*: If the number of distinct measure rearrangements that can alter their initial pair and set distribution is greater than *max_rearrangements*, the function *muMCWtest* performs N = *max_rearrangements* random measure rearrangements to calculate the collection of expected-by-chance muMCW_BIs.
#'
#'    - It calculates the *P~upper~* and *P~lower~* values as the fraction of expected-by-chance muMCW_BIs that are higher or equal to and lower or equal to the observed muMCW_BIs, respectively.
#'
#' @format
#' When executing the *muMCWtest* function, users must provide the path to a local CSV file named *X_muMCWtest_data.csv*, where *X* serves as a user-defined identifier. *X_muMCWtest_data.csv* can be structured in two distinct formats:
#'
#'  - **Vertical layout**: This format allows appending datasets with varying structures, such as different numbers of measure matched-pairs for each appended test. Vertical entry datasets should include the following columns:
#'
#'      - Columns *condition_a* and *condition_b* uniquely identify the two conditions under which matched-paired measures were collected.
#'      - Columns *value_a* and *value_b* contain the actual measures under analysis.
#'      - As many informative columns as needed by users to contextualize the results of each test. The names of the these columns should not contain the terms *condition* or *value*. While these columns are optional when running a single test, at least one column is required when running multiple tests simultaneously. All rows for each individual test must contain the same information in these columns.
#'
#'  - **Horizontal layout**: This format allows appending datasets with similar structures, such as the same number of matched-paired measures for each appended test. Horizontal entry datasets should include the following columns:
#'
#'      - Columns *condition_a* and *condition_b* uniquely identify the two conditions under which matched-paired measures were collected.
#'      - Columns *a.i* and *b.i*, where *i* represents integers to differentiate each specific matched-pairs of measures, contain the actual measures under analysis.
#'      - As many informative columns as needed by users to contextualize the results of each test. The names of these columns should not contain the term *condition* or have the same structure as the *a.i* and *b.i* columns. While these columns are optional when running a single test, at least one column is required when running multiple tests simultaneously.
#'
#' @param path
#' Path for the local CSV file containing the entry dataset formatted for muMCW tests.
#'
#' @param max_rearrangements
#' User-defined maximum number of rearrangements of the dataset used by the function *muMCWtest* to generate a collection of expected-by-chance muMCW_BIs and estimate the statistical significance of observed muMCW_BIs. If the number of distinct dataset rearrangements is less than *max_rearrangements*, *muMCWtest* calculates muMCW_BIs for all possible data rearrangements. If the number of distinct dataset rearrangements is greater than *max_rearrangements*, *muMCWtest* will perform N = *max_rearrangements* random measure rearrangements to calculate the collection of expected-by-chance muMCW_BIs.
#'
#' @return
#' The *muMCWtest* function reports to the console the total number of tests it will execute, and their exact and approximated counts. It also creates a CSV file named *X_muMCWtest_results.csv* where *X* is a user-defined identifier for the entry dataset CSV file. The *X_muMCWtest_results.csv* file contains two rows for each muMCW test, with muMCW_BIs calculated for each possible condition contrast (*e.g.*, *a-b* and *b-a*). The *X_muMCWtest_results.csv* file includes the following columns:
#'
#'    - User-provided informative columns to contextualize the results of each test.
#'    - Columns *condition_a* and *condition_b* indicate the two conditions for which matched-paired measures were provided.
#'    - Column *N* indicates the total number of measure matched-pairs after removing matched-pairs with missing values (NAs).
#'    - Column *test_type* distinguishes between exact and approximated tests.
#'    - Column *BI_type* indicates muMCW_BI.
#'    - Column *condition_contrast* indicates the condition contrast for each row of results.
#'    - Column *observed_BI* contains the value of muMCW_BIs obtained from analyzing the user-provided dataset.
#'    - Column *expected_by_chance_BI_N* indicates the number of data rearrangements used to calculate the expected-by-chance muMCW_BIs. This value corresponds to the lowest number between all possible measure rearrangements and the parameter *max_rearrangements*.
#'    - Columns *pupper* and *plower* represent *P~upper~* and *P~lower~* values, respectively. They denote the fraction of expected-by-chance muMCW_BIs with values higher or equal to and lower or equal to the observed muMCW_BIs, respectively.
#'
#' @import data.table
#' @import pbapply
#' @import Zseq
#'
#' @examples
#' # running muMCWtest with a vertical entry dataset
#' path_v <- system.file("extdata", "muMCWtest_vertical_example_data.csv", package = "MCWtests")
#' muMCWtest_v_results <- muMCWtest(path_v, 200)
#'
#' # running muMCWtest with an horizontal entry dataset
#' path_h <- system.file("extdata", "muMCWtest_horizontal_example_data.csv", package = "MCWtests")
#' muMCWtest_h_results <- muMCWtest(path_h, 200)
#'
#' @export
muMCWtest <- function(path, max_rearrangements) {

  environment(muMCW_approximated_test) <- environment()

  # loading entry dataset using the provided path
  muMCWtest_dataset <- fread(path)
  if (length(grep("value", colnames(muMCWtest_dataset))) == 0) {
    muMCWtest_dataset <- na.omit(dcast(suppressWarnings(melt(muMCWtest_dataset,
                                                             measure.vars = patterns("\\."),
                                                             variable.name = "condition",
                                                             value.name = "value",
                                                             value.factor = FALSE))[
                                                               , `:=`(c("condition", "replicate_pair"), tstrsplit(condition, ".", fixed = TRUE))][
                                                                 , `:=`(replicate_pair = as.numeric(replicate_pair))][
                                                                   , `:=`(condition, paste0("value_", condition))],
                                       ... ~ condition, value.var = "value")[
                                         , `:=`(replicate_pair, NULL)])
  } else if (length(grep("value", colnames(muMCWtest_dataset))) != 0) {
    muMCWtest_dataset <- na.omit(muMCWtest_dataset)
  }

  # determining exact/approximated tests
  muMCWtest_triage <- merge(copy(muMCWtest_dataset),
                            unique(copy(muMCWtest_dataset)[
                              , `:=`(value_a = NULL, value_b = NULL)])[
                                , `:=`(internal_contrast_ID = .I)])[
                                  , `:=`(N = .N), by = internal_contrast_ID][
                                    , `:=`(possible_rearrangements, do.call(rbind, lapply(N, function(x) {
                                      as.numeric(max(Factorial.Double(x+1), gmp = TRUE, odd = TRUE))})))][
                                        , `:=`(test_type, fifelse(is.finite(possible_rearrangements) & possible_rearrangements < max_rearrangements, "exact", "approximated"))][
                                          , `:=`(replicate_pair, seq_len(.N)), by = internal_contrast_ID][
                                            , `:=`(possible_rearrangements, NULL)]

  # returning messages to the console with the total number of tests and the assortment of exact and approximated tests to be run.
  message(paste0("total number of tests: ", NROW(unique(muMCWtest_triage[, internal_contrast_ID])), "\n"),
          paste0("number of exact tests: ", NROW(unique(muMCWtest_triage[test_type == "exact", internal_contrast_ID])), "\n"),
          paste0("number of approximated tests: ", NROW(unique(muMCWtest_triage[test_type == "approximated", internal_contrast_ID]))))

  # running exact and approximated tests
  muMCWtest_results <- setorder(unique(copy(muMCWtest_triage)[
    ,`:=`(grep("replicate_pair|value", colnames(muMCWtest_triage)), NULL)])[
      setcolorder(rbind(if (NROW(unique(muMCWtest_triage[test_type == "exact", internal_contrast_ID])) > 0) {
        copy(muMCWtest_triage)[test_type == "exact", muMCW_exact_test(.SD)]
      }, if (NROW(unique(muMCWtest_triage[test_type == "approximated", internal_contrast_ID])) > 0) {
        copy(muMCWtest_triage)[test_type == "approximated", muMCW_approximated_test(.SD)]
      }), c("internal_contrast_ID", "BI_type", "condition_contrast", "observed_BI", "expected_by_chance_BI_N", "pupper", "plower")),
      on = .(internal_contrast_ID = internal_contrast_ID)][
        , `:=`(internal_contrast_ID, as.numeric(gsub("internal_contrast_", "", internal_contrast_ID)))], internal_contrast_ID)[
          , `:=`(internal_contrast_ID, NULL)]

  # writing results table to the same directory from which the entry dataset was retrieved
  fwrite(muMCWtest_results, file = gsub("muMCWtest_data.csv", "muMCWtest_results.csv", path, fixed = TRUE), quote = FALSE, row.names = FALSE)
  return(muMCWtest_results)
}

utils::globalVariables(codetools::findGlobals(muMCWtest, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(muMCWtest, merge = FALSE)$variables)

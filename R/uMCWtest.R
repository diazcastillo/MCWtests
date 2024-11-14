#' @title
#' Unmatched-Measures Monte Carlo-Wilcoxon (uMCW) Test
#'
#' @description
#' The uMCW test is a statistical tool to assess whether two sets of unmatched measures and their heterogeneity are significantly biased in the same direction. Significantly different data heterogeneities between two conditions could indicate that the measure under analysis is more constrained or more relaxed in one of the conditions, potentially providing insights into the mechanisms underlying the variation of such measure. For instance, uMCW tests can be used to analyze bodyweights or transcript abundances determined for two sets of mice that have been maintained in different conditions.
#’
#' @details
#' The function *uMCWtest* eliminate missing values (NAs) from the dataset before proceeding these steps.
#'
#' - To estimate the bias between the two sets of measures (*e.g.*, *a* and *b*), the function *uMCWtest* performs these tasks:
#'
#'    - It generates all possible disjoint data pairs using measures from both sets.
#'    - For each measure pair, it subtracts the second measure in the pair from the first measure in the pair.
#'    - It ranks the absolute values of all non-zero measure pair differences from lowest to highest. Measure pair differences with a value of 0 are assigned a 0 rank. If multiple measure pair differences have the same absolute value, all tied measure pair differences are assigned the lowest rank possible.
#'    - It assigns each measure pair rank a sign based on the sign of its corresponding measure pair difference.
#'    - It sums the signed ranks for measure pairs formed with measures from the two different sets (*e.g.*, *a-b* and *b-a*).
#'    - For each type of disjoint set measure pairs (*e.g.*, *a-b* and *b-a*), it calculates uMCW_BI by dividing the sum of signed ranks by the maximum number this sum could have if the corresponding measure pairs had the highest possible positive ranks. Consequently, uMCW_BI ranges between 1 when all the values for measures in the first set are higher than all the values from measures in the second set, and -1 when all the values for measures in the first set are lower than all the values from measures in the second set.
#'
#'  - To estimate the bias between the heterogeneity of two sets of measures, the function *uMCWtest* performs these tasks:
#'
#'    - It generates all possible disjoint data pairs within each set, disregarding the order of the paired measures. For instance, the measure pair *a.1-a.2* is considered equivalent to the measure pair *a.2-a.1*, and only the former is retained for the subsequent calculations.
#'    - For each measure pair, it subtracts the second measure from the first measure.
#'    - It ranks all measure pair differences with non-zero values from lowest to highest. Measure pair differences with a value of 0 are assigned a 0 rank. If multiple measure pair differences have the same absolute value, *uMCWtest* assigns all tied measure pair differences the lowest rank possible.
#'    - It sums ranks for measure pairs formed with measures from the same set (*e.g.*, *a-a* and *b-b*).
#'    - For each type of same-set measure pairs (*e.g.*, *a-a* and *b-b*), it divides each sum of signed ranks by the maximum number this sum could have if the corresponding measure pairs had the highest possible ranks.
#'    - It calculates two heterogeneity bias indexes (uMCW_HBIs) by subtracting the normalized sum of signed ranks from the previous step in two possible directions (*e.g.*, *a-b* and *b-a*). Consequently, uMCW_HBI ranges between 1 when at least two measures in the first set have distinct values and all measures in the second set have the same value, and -1 when all measures in the first set have the same value and at least two measures in the second set have distinct values.
#'
#'  - To assess the significance of the uMCW_BIs and uMCW_HBIs obtained with the user-provided data (observed uMCW_BIs and uMCW_HBIs), the function *uMCWtest* performs these tasks:
#'
#'    - It generates a collection of expected-by-chance uMCW_BIs and uMCW_HBIs. These expected values are obtained by rearranging the measures between the two sets multiple times. The user-provided parameter *max_rearrangements* determines the two paths that the function *uMCWtest* can follow to generate the collection of expected-by-chance uMCW_BIs and uMCW_HBIs:
#'
#'      - *uMCW exact testing*: If the number of distinct measure rearrangements that can alter their initial set distribution is less than *max_rearrangements*, the function *uMCWtest* calculates uMCW_BIs and uMCW_HBIs for all possible data rearrangements.
#'      - *uMCW approximated testing*: If the number of distinct measure rearrangements that can alter their initial set distribution is greater than *max_rearrangements*, the function *uMCWtest* will perform N = *max_rearrangements* random measure rearrangements to calculate the collection of expected-by-chance uMCW_BIs and uMCW_HBIs.
#'
#'    - It calculates *P~upper~* and *P~lower~* values as the fraction of expected-by-chance uMCW_BIs and uMCW_HBIs that are higher or equal to and lower or equal to the observed uMCW_BIs and uMCW_HBIs, respectively.
#’
#' @format
#' When executing the *uMCWtest* function, users must provide the path to a local CSV file named *X_uMCWtest_data.csv*, where *X* serves as a user-defined identifier. *X_uMCWtest_data.csv* can be structured in two distinct formats:
#'
#'  - **Vertical layout**: This format allows appending datasets with varying structures, such as different numbers of measures per set or between each appended test. Vertical entry datasets should include the following columns:
#'
#'    - The *condition* column uniquely identifies each of the two measure sets under analysis.
#'    - The *value* column contains the actual measures under analysis.
#'    - As many informative columns as needed by users to contextualize the results of each test. The names of these columns should not include the terms *condition* or *value*. While these columns are optional when running a single test, at least one column is required when running multiple tests simultaneously. All rows for each individual test must contain the same information in these columns.
#'
#'  - **Horizontal layout**: This format allows appending datasets with similar structures, such as the same number of measures collected for each of the two conditions. Horizontal entry datasets should include the following columns:
#'
#'    - Columns *condition_a* and *condition_b* uniquely identify the two measure sets under analysis.
#'    - Columns *a.i* and *b.j*, where *i* and *j* represent integers to differentiate specific measures within each set, contain the actual measures under analysis.
#'    - As many informative columns as needed by users to contextualize the results of each test. The names of these columns should not contain the term *condition* or have the same structure as the *a.i* and *b.j* columns. While these columns are optional when running a single test, at least one column is required when running multiple tests simultaneously.
#'
#' @param path
#' Path for the local CSV file containing the entry dataset formatted for uMCW tests.
#'
#' @param max_rearrangements
#' User-defined maximum number of rearrangements of the dataset used by the function *uMCWtest* to generate a collection of expected-by-chance uMCW_BIs and uMCW_HBIs and estimate the statistical significance of observed uMCW_BIs and uMCW_HBIs. If the number of distinct dataset rearrangements is less than *max_rearrangements*, *uMCWtest* calculates uMCW_BIs and uMCW_HBIs for all possible data rearrangements. If the number of distinct dataset rearrangements is greater than *max_rearrangements*, *uMCWtest* will perform N = *max_rearrangements* random measure rearrangements to calculate the collection of expected-by-chance uMCW_BIs and uMCW_HBIs.
#’
#' @return
#' The *uMCWtest* function reports to the console the total number of tests it will execute, and their exact and approximated counts. It also creates a CSV file named *X_uMCWtest_results.csv*, where *X* is a user-defined identifier for the entry dataset CSV file. The *X_uMCWtest_results.csv* file contains four rows for each uMCWtest, two for uMCW_BIs calculated for each condition contrast (*e.g.*, *a-b* and *b-a*), and two for uMCW_HBIs calculated for each condition contrast. The *X_uMCWtest_results.csv* file includes the following columns:
#'
#'  - User-provided informative columns to contextualize the results of each test.
#'  - Columns *condition_a* and *condition_b* indicate the two measure sets under analysis.
#'  - Columns *N*, *n_a* and *n_b* indicate the total number of measures and the number of measures belonging to each set after removing missing values (NAs).
#'  - Column *test_type* distinguishes between exact and approximated tests.
#'  - Column *BI_type* indicates the bias index type (uMCW_BI and uMCW_HBI) for each row of results.
#'  - Column *condition_contrast* indicates the set contrast (*e.g.*, *a-b* or *b-a*) for each row of results.
#'  - Column *observed_BI* contains the values of uMCW_BIs and uMCW_HBIs obtained from analyzing the user-provided dataset.
#'  - Column *expected_by_chance_BI_N* indicates the number of data rearrangements used to calculate the expected-by-chance uMCW_BIs and uMCW_HBIs. This value corresponds to the lowest number between all possible measure rearrangements and the parameter *max_rearrangements*.
#'  - Columns *pupper* and *plower* represent the *P~upper~* and *P~lower~* values, respectively. They denote the fraction of expected-by-chance uMCW_BIs or uMCW_HBIs with values higher or equal to and lower or equal to the observed uMCW_BIs or uMCW_HBIs, respectively.
#'
#' @import data.table
#' @import pbapply
#'
#' @examples
#' # Executing uMCWtest with an ideal vertical entry dataset
#' path_v <- system.file("extdata", "example_vertical_uMCWtest_data.csv", package = "MCWtests")
#' uMCWtest_v_results <- uMCWtest(path_v, 200)
#'
#' # Executing uMCWtest with an ideal horizontal entry dataset
#' path_h <- system.file("extdata", "example_horizontal_uMCWtest_data.csv", package = "MCWtests")
#' uMCWtest_h_results <- uMCWtest(path_h, 200)
#'
#' @export
uMCWtest <- function(path, max_rearrangements) {

  environment(uMCW_approximated_test) <- environment()

  # loading entry dataset using the provided path
  uMCWtest_dataset <- fread(path)
  if(length(grep("value", colnames(uMCWtest_dataset))) == 0) {
    uMCWtest_dataset <- na.omit(suppressWarnings(melt(fread(path),
                                                      measure.vars = patterns("\\."),
                                                      variable.name = "condition",
                                                      value.name = "value",
                                                      value.factor = FALSE)))[
                                                        , `:=`(c("condition", "replicate"), tstrsplit(condition, ".", fixed = TRUE))][
                                                          , `:=`(replicate, NULL)]
  } else if (length(grep("value", colnames(uMCWtest_dataset))) != 0) {
    uMCWtest_dataset <- na.omit(uMCWtest_dataset)[
      , `:=`(condition_a = unique(condition)[1],
             condition_b = unique(condition)[2]),
      by = setdiff(colnames(uMCWtest_dataset), c("condition", "value"))][
        , `:=`(condition = fifelse(condition == unique(condition_a)[1], "a", "b")),
        by = setdiff(colnames(uMCWtest_dataset), c("condition", "value"))]
  }

  # determining exact/approximated tests
  uMCWtest_triage <- setorder(merge(copy(uMCWtest_dataset),
                                    unique(copy(uMCWtest_dataset)[
                                      , `:=`(condition = NULL, value = NULL)])[
                                        , `:=`(internal_contrast_ID = .I)])[
                                          , `:=`(N = .N), by = internal_contrast_ID][
                                            , `:=`(n = .N), by = .(internal_contrast_ID, condition)][
                                              , `:=`(n_a = unique(n)[1], n_b = fifelse(length(unique(n)) > 1, unique(n)[2], unique(n)[1])), by = internal_contrast_ID][
                                                , `:=`(possible_rearrangements, factorial(N) / (factorial(n_a) * factorial(N - n_a)))][
                                                  , `:=`(test_type, fifelse(is.finite(possible_rearrangements) & possible_rearrangements < max_rearrangements, "exact", "approximated"))][
                                                    , `:=`(replicate = seq_len(.N)), by = internal_contrast_ID][
                                                      , `:=`(n = NULL, possible_rearrangements = NULL)], internal_contrast_ID, condition)

  # returning messages to the console with the total number of tests and the assortment of exact and approximated tests to be run.
  message(paste0("total number of tests: ", NROW(unique(uMCWtest_triage[, internal_contrast_ID])), "\n"),
          paste0("number of exact tests: ", NROW(unique(uMCWtest_triage[test_type == "exact", internal_contrast_ID])), "\n"),
          paste0("number of approximated tests: ", NROW(unique(uMCWtest_triage[test_type == "approximated", internal_contrast_ID]))))

  # running exact and approximated tests
  uMCWtest_results <- setorder(unique(copy(uMCWtest_triage)[
    ,`:=`(grep("^condition$|value|replicate", colnames(uMCWtest_triage)), NULL)])[
      setcolorder(rbind(if (NROW(unique(uMCWtest_triage[test_type == "exact", internal_contrast_ID])) > 0) {
        copy(uMCWtest_triage)[test_type == "exact", uMCW_exact_test(.SD)]
      }, if (NROW(unique(uMCWtest_triage[test_type == "approximated", internal_contrast_ID])) > 0) {
        copy(uMCWtest_triage)[test_type == "approximated", uMCW_approximated_test(.SD)]
      }), c("internal_contrast_ID", "BI_type", "condition_contrast", "observed_BI", "expected_by_chance_BI_N", "pupper", "plower")),
      on = .(internal_contrast_ID = internal_contrast_ID), allow.cartesian = TRUE][
        , `:=`(internal_contrast_ID, as.numeric(gsub("internal_contrast_", "", internal_contrast_ID)))], internal_contrast_ID)[
          , `:=`(internal_contrast_ID, NULL)][
            , `:=`(condition_contrast = fifelse(condition_contrast == "a-b", paste0(condition_a, "-", condition_b), paste0(condition_b, "-", condition_a)))]

  # writing results table to the same directory from which the entry dataset was retrieved
  fwrite(uMCWtest_results, file = gsub("uMCWtest_data.csv", "uMCWtest_results.csv", path, fixed = TRUE), quote = FALSE, row.names = FALSE)
  return(uMCWtest_results)
}

utils::globalVariables(codetools::findGlobals(uMCWtest, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(uMCWtest, merge = FALSE)$variables)

#' @title
#' Matched-Measures Bivariate Monte Carlo-Wilcoxon (mbMCW) Test
#'
#' @description
#' The mbMCW test is a statistical tool to assess whether two sets of inherently matched-paired measures are significantly differentially biased in the same direction. For instance, mbMCW tests can be used to analyze bodyweights or transcript abundances determined at two different timepoints for two sets of mice that have been exposed to different conditions.
#'
#' @details
#' The function *mbMCWtest* eliminates any matched-paired measures with at least one missing value (NA) before proceeding with the following steps.
#'
#'  - To estimate the differential bias between the two sets of matched-paired measures in the dataset, the function *mbMCWtest* perfoms the following tasks:
#'
#'      - For each matched-paired measure, it subtracts the values for the two possible matched condition contrasts (*e.g.*, *a-b* and *b-a*).
#'      - For each matched condition contrast, it ranks the absolute values of non-zero differences from lowest to highest. Measure pair differences with a value of 0 are assigned a 0 rank. If multiple measure pair differences have the same absolute value, all tied measure pair differences are assigned the lowest rank possible.
#'      - It assigns each measure pair rank a sign based on the sign of its corresponding measure pair difference.
#'      - For each set of matched-paired measures (*e.g.*, *x* and *y*), it sums the signed ranks for each matched condition contrast (*e.g.*, *a-b* and *b-a*).
#'      - For each set of matched-paired measures (*e.g.*, *x* and *y*) and each matched condition contrast (*e.g.*, *a-b* and *b-a*), it calculates one mbMCW_BI. This value is obtained by dividing each sum of signed ranks by the maximum number this sum could have if the corresponding measure pairs had the highest possible positive ranks. Consequently, mbMCW_BI ranges between 1 when all the values for matched-pair measure differences in the set under analysis have the highest positive values, and -1 when all the values for matched-pair measure differences in the set under analysis have the lowest negative values.
#'
#'  - To assess the significance of the mbMCW_BIs obtained from the user-provided dataset (observed mbMCW_BIs), the function *mbMCWtest* perfoms the following tasks:
#'
#'      - It generates a collection of expected-by-chance mbMCW_BIs. These expected values are obtained by rearranging the matched-pair measures between the two sets multiple times. The user-provided parameter *max_rearrangements* determines the two paths the function *mbMCWtest* can follow to generate the collection of expected-by-chance mbMCW_BIs:
#'
#'        - *mbMCW exact testing*: If the number of distinct matched-paired measure rearrangements that can alter their initial set distribution is less than *max_rearrangements*, the function *mbMCWtest* calculates mbMCW_BIs for all possible data rearrangements.
#'        - *mbMCW approximated testing*: If the number of distinct matched-paired measure rearrangements that can alter their initial set distribution is greater than *max_rearrangements*, the function *mbMCWtest* performs N = *max_rearrangements* random measure rearrangements to calculate the collection of expected-by-chance mbMCW_BIs.
#'
#'      - It calculates *P~upper~* and *P~lower~* values as the fraction of expected-by-chance mbMCW_BIs that are higher or equal to and lower or equal to the observed mbMCW_BIs, respectively.
#'
#' @format
#' When executing the *mbMCWtest* function, users must provide the path to a local CSV file named *X_mbMCWtest_data.csv*, where *X* serves as a user-defined identifier. *X_mbMCWtest_data.csv* can be structured in two distinct formats:
#'
#'  - **Vertical layout:** This format allows appending datasets with varying structures, such as different numbers of matched-pairs per set or between each appended test. Vertical entry datasets should include the following columns:
#'
#'      - Columns *matched_condition_a* and *matched_condition_b* uniquely identify the two conditions under which matched-paired measure were collected.
#'      - Column *unmatched_condition* uniquely identifies the two sets of matched-paired measures under analysis.
#'      - Columns *value_a* and *value_b* contain the actual measures under analysis.
#'      - As many informative columns as needed by users to contextualize the results of each test. The names of these columns should not contain the terms *condition* or *value*. While these columns are optional when running a single test, at least one column is required when running multiple tests simultaneously. All rows for each individual test must contain the same information in these columns.
#'
#'  - **Horizontal layout:** This format allows appending datasets with similar structures, such as the same number of matched-paired measures collected for two conditions. Horizontal entry datasets should include the following columns:
#'
#'      - Columns *matched_condition_a* and *matched_condition_b* uniquely identify the two conditions under which matched-paired measure were collected.
#'      - Columns *unmatched_condition_x* and *unmatched_condition_y* uniquely identify the two different sets of matched-paired measures under analysis.
#'      - Columns *x.a.i*, *y.a.i*, *x.b.i* and *y.b.i*, where *i* represents integers to differentiate each specific matched-pair of measures, contain the actual measures under analysis.
#'      - As many informative columns as needed by users to contextualize the results of each test. The name of these columns should not contain the term *condition* or have the same structure as the *x.a.i*, *y.a.i*, *x.b.i* and *y.b.i* columns. While these columns are optional when running a single test, at least one column is required when running multiple tests simultaneously.
#'
#' @param path
#' Path for the local CSV file containing the entry dataset formatted for mbMCW tests.
#'
#' @param max_rearrangements
#' User-defined maximum number of rearrangements of the dataset used by the function *mbMCWtest* to generate a collection of expected-by-chance mbMCW_BIs and estimate the statistical significance of observed mbMCW_BIs. If the number of distinct dataset rearrangements is less than *max_rearrangements*, *mbMCWtest* calculates mbMCW_BIs for all possible data rearrangements. If the number of distinct dataset rearrangements is greater than *max_rearrangements*, *mbMCWtest* will perform N = *max_rearrangements* random measure rearrangements to calculate the collection of expected-by-chance mbMCW_BIs.
#'
#' @return
#' The *mbMCWtest* function reports to the console the total number of tests it will execute, and their exact and approximated counts. It also creates a CSV file named *X_mbMCWtest_results.csv*, where *X* is a user-defined identifier for the entry dataset CSV file. The *X_mbMCWtest_results.csv* file contains four rows for each mbMCWtest, with mbMCW_BIs calculated for each possible contrast between matched and unmatched measures (*e.g.*, *a-b*, *b-a*, *x-y* and *y-x*). The *X_mbMCWtest_results.csv* file includes the following columns:
#'
#'    - User-provided informative columns to contextualize the results of each test.
#'    - Columns *matched_condition_a* and *matched_condition_b* indicate the conditions for which matched-paired measures were provided.
#'    - Columns *unmatched_condition_x* and *unmatched_condition_y* indicate the two sets of matched-pairs measures.
#'    - Columns *N*, *N_x* and *N_y* indicate the total number of matched-paired measures, and their distribution between the two unmatched sets after removing any matched-pair with missing values (NAs).
#'    - Column *test_type* distinguishes between exact and approximated tests.
#'    - Column *BI_type* indicates mbMCW_BI.
#'    - Column *matched_condition_contrast* and *unmatched_condition_contrast* indicate the matched and unmatched condition contrast for each row of results.
#'    - Column *observed_BI* contains the value of mbMCW_BIs obtained from analyzing the user-provided dataset.
#'    - Column *expected_by_chance_BI_N* indicates the number of data rearrangements used to calculate the expected-by-chance mbMCW_BIs. This value corresponds to the lowest number between all possible measure rearrangements and the parameter *max_rearrangements*.
#'    - Columns *pupper* and *plower* represent the *P~upper~* and *P~lower~* values, respectively. They denote the fraction of expected-by-chance mbMCW_BIs with values higher or equal to and lower or equal to the observed mbMCW_BIs, respectively.
#'
#' @import data.table
#' @import pbapply
#'
#' @examples
#'
#' test_temp <- tempdir()
#' extdata_v <- system.file("extdata", "example_vertical_mbMCWtest_data.csv", package = "MCWtests")
#' file.copy(extdata_v, test_temp)
#' extdata_h <- system.file("extdata", "example_horizontal_mbMCWtest_data.csv", package = "MCWtests")
#' file.copy(extdata_h, test_temp)
#' # running mbMCWtest with an ideal vertical entry dataset
#' path_v <- file.path(test_temp, "example_vertical_mbMCWtest_data.csv")
#' mbMCWtest_vertical_results <- mbMCWtest(path_v, 10)
#' print(mbMCWtest_vertical_results)
#' # running mbMCWtest with an ideal horizontal entry dataset
#' path_h <- file.path(test_temp, "example_horizontal_mbMCWtest_data.csv")
#' mbMCWtest_horizontal_results <- mbMCWtest(path_h, 10)
#' print(mbMCWtest_horizontal_results)
#'
#' rm(test_temp)
#'
#' @export
mbMCWtest <- function(path, max_rearrangements) {

  environment(mbMCW_approximated_test) <- environment()

  # loading entry dataset using the provided path
  mbMCWtest_dataset <- fread(path)
  if (length(grep("value", colnames(mbMCWtest_dataset))) == 0) {
    mbMCWtest_dataset <- na.omit(dcast(suppressWarnings(melt(mbMCWtest_dataset,
                                                             measure.vars = patterns("\\."),
                                                             variable.name = "condition",
                                                             value.name = "value",
                                                             value.factor = FALSE))[
                                                               value != "NA", ][
                                                                 , `:=`(c("unmatched_condition", "matched_condition", "replicate_pair"), tstrsplit(condition, ".", fixed = TRUE))][
                                                                   , `:=`(matched_condition = paste0("value_", matched_condition),
                                                                          condition = NULL)],
                                       ... ~ matched_condition, value.var = "value"))[
                                         , `:=`(replicate_pair, NULL)]
  } else if (length(grep("value", colnames(mbMCWtest_dataset))) != 0) {
    mbMCWtest_dataset <- na.omit(mbMCWtest_dataset)[
      , `:=`(unmatched_condition_x = unique(unmatched_condition)[1],
             unmatched_condition_y = unique(unmatched_condition)[2]),
      by = setdiff(colnames(mbMCWtest_dataset), c("replicate_pair", "unmatched_condition", "value_a", "value_b"))][
        , `:=`(unmatched_condition = fifelse(unmatched_condition == unmatched_condition_x, "x", "y"))]
  }

  # determining exact/approximated tests
  mbMCWtest_triage <- merge(copy(mbMCWtest_dataset),
                            unique(copy(mbMCWtest_dataset)[
                              , `:=`(unmatched_condition = NULL, value_a = NULL, value_b = NULL)])[
                                , `:=`(internal_contrast_ID = .I)])[
                                  , `:=`(N = .N), by = internal_contrast_ID][
                                    , `:=`(n = .N), by = .(internal_contrast_ID, unmatched_condition)][
                                      , `:=`(N_x = unique(n)[1], N_y = fifelse(length(unique(n)) > 1, unique(n)[2], unique(n)[1])), by = internal_contrast_ID][
                                        , `:=`(possible_rearrangements, factorial(N) / (factorial(n) * factorial(N - n)))][
                                          , `:=`(test_type, fifelse(is.finite(possible_rearrangements) & possible_rearrangements < max_rearrangements, "exact", "approximated"))][
                                            , `:=`(n = NULL, possible_rearrangements = NULL)][
                                              , `:=`(`a-b_condition_contrast` = "a-b",
                                                     `b-a_condition_contrast` = "b-a",
                                                     abs_diff = abs(value_a - value_b))][
                                                       abs_diff == 0, `:=`(abs_diff = NA)][
                                                         , `:=`(abs_rank = rank(abs_diff, na.last = "keep", ties.method = "min")), by = internal_contrast_ID][
                                                           is.na(abs_rank), `:=`(abs_rank = 0)][
                                                             , `:=`(`a-b_signed_rank` = abs_rank * sign(value_a - value_b),
                                                                    `b-a_signed_rank` = abs_rank * sign(value_b - value_a))][
                                                                      , `:=`(replicate_pair, seq_len(.N)), by = internal_contrast_ID]

  # returning messages to the console with the total number of tests and the assortment of exact and approximated tests to be run.
  message(paste0("total number of tests: ", NROW(unique(mbMCWtest_triage[, internal_contrast_ID])), "\n"),
          paste0("number of exact tests: ", NROW(unique(mbMCWtest_triage[test_type == "exact", internal_contrast_ID])), "\n"),
          paste0("number of approximated tests: ", NROW(unique(mbMCWtest_triage[test_type == "approximated", internal_contrast_ID]))))

  # running exact and approximated tests
  mbMCWtest_results <- setorder(unique(copy(mbMCWtest_triage)[
    ,`:=`(grep("^unmatched_condition$|replicate_pair|value_|a-b|b-a|abs_", colnames(mbMCWtest_triage)), NULL)])[
      setcolorder(rbind(if (NROW(unique(mbMCWtest_triage[test_type == "exact", internal_contrast_ID])) > 0) {
        copy(mbMCWtest_triage)[test_type == "exact", mbMCW_exact_test(.SD)]
      }, if (NROW(unique(mbMCWtest_triage[test_type == "approximated", internal_contrast_ID])) > 0) {
        copy(mbMCWtest_triage)[test_type == "approximated", mbMCW_approximated_test(.SD)]
      }), c("internal_contrast_ID", "BI_type", "matched_condition_contrast", "unmatched_condition_contrast", "observed_BI", "expected_by_chance_BI_N", "pupper", "plower")),
      on = .(internal_contrast_ID), allow.cartesian = TRUE], internal_contrast_ID)[
        , `:=`(internal_contrast_ID, NULL)][
          , `:=`(matched_condition_contrast = fifelse(matched_condition_contrast == "a-b", paste0(matched_condition_a, "-", matched_condition_b), paste0(matched_condition_b, "-", matched_condition_a)),
                 unmatched_condition_contrast = fifelse(unmatched_condition_contrast == "x-y", paste0(unmatched_condition_x, "-", unmatched_condition_y), paste0(unmatched_condition_y, "-", unmatched_condition_x)))]

  # writing results table to the same directory from which the entry dataset was retrieved
  fwrite(mbMCWtest_results, file = gsub("mbMCWtest_data.csv", "mbMCWtest_results.csv", path, fixed = TRUE), quote = FALSE, row.names = FALSE)
  return(mbMCWtest_results)
}

utils::globalVariables(codetools::findGlobals(mbMCWtest, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(mbMCWtest, merge = FALSE)$variables)

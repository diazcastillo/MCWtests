#' @title
#' Bias-measures Monte Carlo-Wilcoxon (bMCW) Test
#'
#' @description
#' The bMCW test is a statistical tool to assess whether a set of measures of bias for a quantitative trait between two conditions or a subset of these bias measures are themselves significantly biased in the same direction. For instance, bMCW tests can be used to analyze bias indexes obtained using other MCW tests or fold change for transcript abundances spanning the entire transcriptome or only for genes located in specific genomic regions from two sets of mice exposed to different conditions.
#'
#' @details
#' The function *bMCWtest* eliminates missing values (NAs) from the dataset before proceeding with the following steps.
#'
#'  - To estimate the bias for all bias measures in the entire dataset or a subset of them, the function *bMCWtest* performs the following tasks:
#'
#'      - It ranks all bias measures with non-zero values from lowest to highest. Bias measures with a value of 0 are assigned a 0 rank. If multiple bias measures have the same absolute value, all tied bias measures are asssigned the lowest rank possible.
#'      - It assigns each rank a sign based on the sign of its corresponding bias measure.
#'      - It calculates a whole-set bias index (bMCW_wBI) by summing the signed ranks for all elements in the dataset and dividing it by the maximum number that sum could have if all bias measures were positive. Consequently, bMCW_wBI ranges between 1 when all bias measures are positive, and -1 when all bias measures are negative.
#'      - It calculates a subset bias index (bMCW_sBI) for each subset of elements under analysis by summing the signed ranks for the elements in the subset and dividing it by the maximum number that sum could have if the elements in the subset had the highest possible positive bias measures. Consequently, bMCW_sBI ranges between 1 when the bias measures for the subset in question have the highest positive bias measures in the entire dataset, and -1 when the bias measures for the subset in question have the lowest negative bias measures in the entire dataset.
#'
#'  - To assess the significance of the bMCW-wBIs and bMCW-sBIs obtained from the user-provided dataset (observed bMCW-wBIs and bMCW-sBIs), the function *bMCWtest* performs the following tasks,
#'
#'      - It generates a collection of expected-by-chance bMCW_wBIs by rearranging the signs of all signed ranks multiple times. The function *bMCWtest* also generates a collection of expected-by-chance bMCW_sBIs by rearranging the subset of elements multiple times. The user-provided parameter *max_rearrangements* determines the two paths that the function *bMCWtest* can follow to generate the collection of expected-by-chance bMCW_wBIs and bMCW_sBIs:
#'
#'          - *bMCW exact testing*: If the number of distinct bias measure rearrangements that can alter their initial sign distribution or subset distribution is less than *max_rearrangements*, the function *bMCWtest* calculates bMCW_wBIs or bMCW_sBIs for all possible data rearrangements.
#'          - *bMCW approximated testing*: If the number of distinct bias measure rearrangements that can alter their initial sign distribution or subset distribution is greater than *max_rearrangements*, the function *bMCWtest* performs N = *max_rearrangements* random measure rearrangements to calculate the collection of expected-by-chance bMCW_wBIs or bMCW_sBIs.
#'
#'      - It calculates *P~upper~* and *P~lower~* values, as the fraction of expected-by-chance bMCW-wBIs and bMCW-sBIs that are higher or equal to and lower or equal to the observed bMCW-wBIs and bMCW-sBIs, respectively.
#’
#' @format
#' When executing the *bMCWtest* function, users must provide the path to a local CSV file named *X_bMCWtest_data.csv*, where *X* serves as a user-defined identifier. *X_bMCWtest_data.csv* should include the following columns:
#'
#'  - Column *bias_value* contains the value of the bias measure under analysis.
#'  - Columns *subset_x*, where *x* represents the specific type of subset for each column, such as "chr" for chromosomes or "GO" for Gene Ontology. These columns are required if users intend to assess whether bias measures for certain subsets of elements in the dataset are significantly biased in the same direction. Columns *subset_x* can indicate whether an element belongs to a subset using either "YES" and "NO", or specific subset names like "chr1" or "chrX", or a combination of both, such as "chr1", "chrX" and "NO". The function *bMCW test* will transform the dataset to conduct independent analysis of each subset of elements marked as "YES" or with a specific subset name in each *subset_x* column.
#'  - As many informative columns as needed by users to contextualize the results of each test. The names of these columns should not contain the terms *bias_value* or *subset*. While these columns are optional when running a single test, at least one column is required when running multiple tests simultaneously. All rows for each individual test must contain the same information in these columns.
#'  - Users can specify columns with information relevant about each element or row using the column name structure *element_x*, where *x* indicates the specific information in each column (see example). However, *element_x* columns are not essential for bMCW testing and will not be included in the results file.
#'
#' @param path
#' Path for the local CSV file containing the entry dataset formatted for bMCW tests.
#'
#' @param max_rearrangements
#' User-defined maximum number of rearrangements of the dataset used by the function *bMCWtest* to generate a collection of expected-by-chance bMCW_wBIs and bMCW_sBIs and estimate the statistical significance of observed bMCW_BIs. If the number of distinct dataset rearrangements is less than *max_rearrangements*, *bMCWtest* calculates bMCW_wBIs and bMCW_sBIs for all possible data rearrangements. If the number of distinct dataset rearrangements is greater than *max_rearrangements*, *bMCWtest* will perform N = *max_rearrangements* random measure rearrangements to calculate the collection of expected-by-chance bMCW_wBIs and bMCW_sBIs.
#'
#' @return
#' The *bMCWtest* function reports to the console the total number of tests it will execute, and their exact and approximated counts. It also creates a CSV file named *X_bMCWtest_results.csv*, where *X* is a user-defined identifier for the entry dataset CSV file. The *X_bMCWtest_results.csv* file contains one row for each bMCWtest to indicate the results of whole-set bMCW testing, and as many rows as necessary to indicate the results of subset bMCW testing. Rows for whole-set analyses will be at the top of *X_bMCWtest_results.csv* file. The *X_bMCWtest_results.csv* file includes the following columns:
#'
#'  - User-provided informative columns to contextualize the results of each test.
#'  - Column *subset_type* indicates whether the results in each row corresponds to *whole-set* tests or specific *subset* tests, such as "chr" for chromosomes or "GO" for Gene Ontology terms.
#'  - Column *tested_subset* indicates the name of the subset under analysis. For *whole-set* tests, the *tested_subset* column indicates "none". For subset tests, the *tested_subset* column indicates "YES" or the specific name of the subset under analysis, such as "chr1" or "chrX".
#'  - Columns *N* and *n* indicate the total number of elements in the whole set and those associated with the subset under analysis, respectively, after removing missing values (NAs). For *whole-set* tests, columns *N* and *n* have the same value.
#'  - Column *test_type* distinguishes between exact and approximated tests.
#'  - Column *BI_type* indicates whether results correspond to *whole-set* tests (bMCW_wBI) or to *subset* tests (bMCW_sBIs).
#'  - Column *observed_BI* contains the value of bMCW_BIs obtained from analyzing the user-provided dataset.
#'  - Column *expected_by_chance_BI_N* indicates the number of data rearrangements used to calculate the expected-by-chance bMCW_wBIs and bMCW_sBIs. This value corresponds to the lowest number between all possible measure rearrangements and the parameter *max_rearrangements*.
#'  - Columns *pupper* and *plower* represent the *P~upper~* and *P~lower~* values, respectively. They denote the fraction of expected-by-chance bMCW_wBIs or bMCW_sBIs with values higher or equal to and lower or equal to the observed bMCW_wBIs or bMCW_sBIs, respectively.
#'
#' @import data.table
#' @import pbapply
#'
#' @examples
#' # running bMCWtest with an ideal entry dataset
#' path <- system.file("extdata", "example_bMCWtest_data.csv", package = "MCWtests")
#' bMCWtest_results <- bMCWtest(path, 200)
#'
#' @export
bMCWtest <- function(path, max_rearrangements) {

  environment(bMCW_approximated_test) <- environment()

  # loading entry dataset using the provided path and determining exact/approximated tests
  bMCWtest_dataset <- na.omit(copy(fread(path)))
  if (length(grep("subset_", colnames(bMCWtest_dataset))) == 0) {
    bMCWtest_triage <- merge(copy(bMCWtest_dataset),unique(copy(bMCWtest_dataset)[
      , `:=`(grep("element_|bias_value", colnames(bMCWtest_dataset)), NULL)])[
        , `:=`(subset_type = "wholeset", tested_subset = "none", internal_contrast_ID = .I)])[
          , `:=`(N = .N), by = internal_contrast_ID][
            , `:=`(test_type = ifelse(max_rearrangements < factorial(N), "approximated", "exact"))][
              , `:=`(internal_element_ID = seq_len(.N),
                     element_rank = paste0("r", seq_len(.N)),
                     element_sign = paste0("s", seq_len(.N)))
              , by = internal_contrast_ID][
                , `:=`(abs_bias_value = abs(bias_value),
                       sign_bias_value = sign(bias_value))][
                         abs_bias_value == 0, `:=`(abs_bias_value, NA)][
                           , `:=`(rank_abs_bias_value = as.numeric(rank(abs_bias_value, na.last = "keep", ties.method = "min"))), by = internal_contrast_ID][
                             is.na(rank_abs_bias_value), `:=`(rank_abs_bias_value, 0)]
  } else if (length(grep("subset_", colnames(bMCWtest_dataset))) != 0) {
    bMCWtest_triage <- merge(rbind(unique(copy(bMCWtest_dataset)[
      , `:=`(grep("element_|subset_|bias_value", colnames(bMCWtest_dataset)), NULL)])[
        , `:=`(subset_type = "wholeset", tested_subset = "none")],
      unique(melt(copy(bMCWtest_dataset)[
        , `:=`(grep("element_|bias_value", colnames(bMCWtest_dataset)), NULL)],
        measure.vars = patterns("subset_"),
        variable.name = "subset_type",
        value.name = "tested_subset"))[
          tested_subset != "NO",])[
            , `:=`(internal_contrast_ID, .I)],
      rbind(copy(bMCWtest_dataset)[
        , `:=`(grep("subset", colnames(bMCWtest_dataset)), NULL)][
          , `:=`(subset_type = "wholeset", subset = "none")],
        melt(copy(bMCWtest_dataset),
             measure.vars = patterns("subset_"),
             variable.name = "subset_type",
             value.name = "subset")),
      allow.cartesian = TRUE)[
        tested_subset != subset, `:=`(subset = "NO")][
          , `:=`(N = .N), by = internal_contrast_ID][
            , `:=`(n = .N), by = .(internal_contrast_ID, subset)][
              , `:=`(possible_rearrangements = fcase(subset_type == "wholeset", factorial(N),
                                                     subset_type != "wholeset", factorial(N) / (factorial(n) * factorial(N - n))))][
                                                       , `:=`(test_type = fifelse(is.finite(possible_rearrangements) & possible_rearrangements < max_rearrangements, "exact", "approximated"))][
                                                         , `:=`(possible_rearrangements, NULL)][
                                                           , `:=`(internal_element_ID = seq_len(.N),
                                                                  element_rank = paste0("r", seq_len(.N)),
                                                                  element_sign = paste0("s", seq_len(.N)))
                                                           , by = internal_contrast_ID][
                                                             , `:=`(abs_bias_value = abs(bias_value),
                                                                    sign_bias_value = sign(bias_value))][
                                                                      abs_bias_value == 0, `:=`(abs_bias_value, NA)][
                                                                        , `:=`(rank_abs_bias_value = as.numeric(rank(abs_bias_value, na.last = "keep", ties.method = "min"))), by = internal_contrast_ID][
                                                                          is.na(rank_abs_bias_value), `:=`(rank_abs_bias_value, 0)]}

  # returning messages to the console with the total number of tests and the assortment of exact and approximated tests to be run.
  message(paste0("total number of tests: ", NROW(unique(bMCWtest_triage[, .(internal_contrast_ID, subset_type)])), "\n"),
          paste0("number of wholeset exact tests: ", NROW(unique(bMCWtest_triage[test_type == "exact" & subset_type == "wholeset", .(internal_contrast_ID, subset_type)])), "\n"),
          paste0("number of subset exact tests: ", NROW(unique(bMCWtest_triage[test_type == "exact" & subset_type != "wholeset", .(internal_contrast_ID, subset_type)])), "\n"),
          paste0("number of wholeset approximated tests: ", NROW(unique(bMCWtest_triage[test_type == "approximated" & subset_type == "wholeset", .(internal_contrast_ID, subset_type)])), "\n"),
          paste0("number of subset approximated tests: ", NROW(unique(bMCWtest_triage[test_type == "approximated" & subset_type != "wholeset", .(internal_contrast_ID, subset_type)]))))

  # running exact and approximated tests
  if(length(grep("^subset$", colnames(bMCWtest_triage))) == 0) {
    bMCWtest_results <- unique(copy(bMCWtest_triage)[,`:=`(grep("element_|^subset$|bias_value", colnames(bMCWtest_triage)), NULL)])
  } else if (length(grep("^subset$", colnames(bMCWtest_triage))) != 0) {
    bMCWtest_results <- unique(copy(bMCWtest_triage)[
      subset != "NO", ][
        ,`:=`(grep("element_|^subset$|bias_value", colnames(bMCWtest_triage)), NULL)])
  }
  bMCWtest_results <- setorder(bMCWtest_results[
    setcolorder(rbind(if (NROW(unique(bMCWtest_triage[test_type == "exact", internal_contrast_ID])) > 0) {
      copy(bMCWtest_triage)[test_type == "exact", bMCW_exact_test(.SD)]
    }, if (NROW(unique(bMCWtest_triage[test_type == "approximated", internal_contrast_ID])) > 0) {
      copy(bMCWtest_triage)[test_type == "approximated", bMCW_approximated_test(.SD)]
    }), c("internal_contrast_ID", "BI_type", "observed_BI", "expected_by_chance_BI_N", "pupper", "plower")),
    on = .(internal_contrast_ID = internal_contrast_ID)][
      , `:=`(internal_contrast_ID, as.numeric(gsub("internal_contrast_", "", internal_contrast_ID)))], internal_contrast_ID)[
        , `:=`(internal_contrast_ID, NULL)]

  # writing results table to the same directory from which the entry dataset was retrieved
  fwrite(bMCWtest_results, file = gsub("bMCWtest_data.csv", "bMCWtest_results.csv", path, fixed = TRUE), quote = FALSE, row.names = FALSE)
  return(bMCWtest_results)
}

utils::globalVariables(codetools::findGlobals(bMCWtest, merge = FALSE)$functions)
utils::globalVariables(codetools::findGlobals(bMCWtest, merge = FALSE)$variables)

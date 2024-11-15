% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbMCWtest.R
\name{mbMCWtest}
\alias{mbMCWtest}
\title{Matched-Measures Bivariate Monte Carlo-Wilcoxon (mbMCW) Test}
\format{
When executing the \emph{mbMCWtest} function, users must provide the path to a local CSV file named \emph{X_mbMCWtest_data.csv}, where \emph{X} serves as a user-defined identifier. \emph{X_mbMCWtest_data.csv} can be structured in two distinct formats:
\itemize{
\item \strong{Vertical layout:} This format allows appending datasets with varying structures, such as different numbers of matched-pairs per set or between each appended test. Vertical entry datasets should include the following columns:
\itemize{
\item Columns \emph{matched_condition_a} and \emph{matched_condition_b} uniquely identify the two conditions under which matched-paired measure were collected.
\item Column \emph{unmatched_condition} uniquely identifies the two sets of matched-paired measures under analysis.
\item Columns \emph{value_a} and \emph{value_b} contain the actual measures under analysis.
\item As many informative columns as needed by users to contextualize the results of each test. The names of these columns should not contain the terms \emph{condition} or \emph{value}. While these columns are optional when running a single test, at least one column is required when running multiple tests simultaneously. All rows for each individual test must contain the same information in these columns.
}
\item \strong{Horizontal layout:} This format allows appending datasets with similar structures, such as the same number of matched-paired measures collected for two conditions. Horizontal entry datasets should include the following columns:
\itemize{
\item Columns \emph{matched_condition_a} and \emph{matched_condition_b} uniquely identify the two conditions under which matched-paired measure were collected.
\item Columns \emph{unmatched_condition_x} and \emph{unmatched_condition_y} uniquely identify the two different sets of matched-paired measures under analysis.
\item Columns \emph{x.a.i}, \emph{y.a.i}, \emph{x.b.i} and \emph{y.b.i}, where \emph{i} represents integers to differentiate each specific matched-pair of measures, contain the actual measures under analysis.
\item As many informative columns as needed by users to contextualize the results of each test. The name of these columns should not contain the term \emph{condition} or have the same structure as the \emph{x.a.i}, \emph{y.a.i}, \emph{x.b.i} and \emph{y.b.i} columns. While these columns are optional when running a single test, at least one column is required when running multiple tests simultaneously.
}
}
}
\usage{
mbMCWtest(path, max_rearrangements)
}
\arguments{
\item{path}{Path for the local CSV file containing the entry dataset formatted for mbMCW tests.}

\item{max_rearrangements}{User-defined maximum number of rearrangements of the dataset used by the function \emph{mbMCWtest} to generate a collection of expected-by-chance mbMCW_BIs and estimate the statistical significance of observed mbMCW_BIs. If the number of distinct dataset rearrangements is less than \emph{max_rearrangements}, \emph{mbMCWtest} calculates mbMCW_BIs for all possible data rearrangements. If the number of distinct dataset rearrangements is greater than \emph{max_rearrangements}, \emph{mbMCWtest} will perform N = \emph{max_rearrangements} random measure rearrangements to calculate the collection of expected-by-chance mbMCW_BIs.}
}
\value{
The \emph{mbMCWtest} function reports to the console the total number of tests it will execute, and their exact and approximated counts. It also creates a CSV file named \emph{X_mbMCWtest_results.csv}, where \emph{X} is a user-defined identifier for the entry dataset CSV file. The \emph{X_mbMCWtest_results.csv} file contains four rows for each mbMCWtest, with mbMCW_BIs calculated for each possible contrast between matched and unmatched measures (\emph{e.g.}, \emph{a-b}, \emph{b-a}, \emph{x-y} and \emph{y-x}). The \emph{X_mbMCWtest_results.csv} file includes the following columns:
\itemize{
\item User-provided informative columns to contextualize the results of each test.
\item Columns \emph{matched_condition_a} and \emph{matched_condition_b} indicate the conditions for which matched-paired measures were provided.
\item Columns \emph{unmatched_condition_x} and \emph{unmatched_condition_y} indicate the two sets of matched-pairs measures.
\item Columns \emph{N}, \emph{N_x} and \emph{N_y} indicate the total number of matched-paired measures, and their distribution between the two unmatched sets after removing any matched-pair with missing values (NAs).
\item Column \emph{test_type} distinguishes between exact and approximated tests.
\item Column \emph{BI_type} indicates mbMCW_BI.
\item Column \emph{matched_condition_contrast} and \emph{unmatched_condition_contrast} indicate the matched and unmatched condition contrast for each row of results.
\item Column \emph{observed_BI} contains the value of mbMCW_BIs obtained from analyzing the user-provided dataset.
\item Column \emph{expected_by_chance_BI_N} indicates the number of data rearrangements used to calculate the expected-by-chance mbMCW_BIs. This value corresponds to the lowest number between all possible measure rearrangements and the parameter \emph{max_rearrangements}.
\item Columns \emph{pupper} and \emph{plower} represent the \emph{P~upper~} and \emph{P~lower~} values, respectively. They denote the fraction of expected-by-chance mbMCW_BIs with values higher or equal to and lower or equal to the observed mbMCW_BIs, respectively.
}
}
\description{
The mbMCW test is a statistical tool to assess whether two sets of inherently matched-paired measures are significantly differentially biased in the same direction. For instance, mbMCW tests can be used to analyze bodyweights or transcript abundances determined at two different timepoints for two sets of mice that have been exposed to different conditions.
}
\details{
The function \emph{mbMCWtest} eliminates any matched-paired measures with at least one missing value (NA) before proceeding with the following steps.
\itemize{
\item To estimate the differential bias between the two sets of matched-paired measures in the dataset, the function \emph{mbMCWtest} perfoms the following tasks:
\itemize{
\item For each matched-paired measure, it subtracts the values for the two possible matched condition contrasts (\emph{e.g.}, \emph{a-b} and \emph{b-a}).
\item For each matched condition contrast, it ranks the absolute values of non-zero differences from lowest to highest. Measure pair differences with a value of 0 are assigned a 0 rank. If multiple measure pair differences have the same absolute value, all tied measure pair differences are assigned the lowest rank possible.
\item It assigns each measure pair rank a sign based on the sign of its corresponding measure pair difference.
\item For each set of matched-paired measures (\emph{e.g.}, \emph{x} and \emph{y}), it sums the signed ranks for each matched condition contrast (\emph{e.g.}, \emph{a-b} and \emph{b-a}).
\item For each set of matched-paired measures (\emph{e.g.}, \emph{x} and \emph{y}) and each matched condition contrast (\emph{e.g.}, \emph{a-b} and \emph{b-a}), it calculates one mbMCW_BI. This value is obtained by dividing each sum of signed ranks by the maximum number this sum could have if the corresponding measure pairs had the highest possible positive ranks. Consequently, mbMCW_BI ranges between 1 when all the values for matched-pair measure differences in the set under analysis have the highest positive values, and -1 when all the values for matched-pair measure differences in the set under analysis have the lowest negative values.
}
\item To assess the significance of the mbMCW_BIs obtained from the user-provided dataset (observed mbMCW_BIs), the function \emph{mbMCWtest} perfoms the following tasks:
\itemize{
\item It generates a collection of expected-by-chance mbMCW_BIs. These expected values are obtained by rearranging the matched-pair measures between the two sets multiple times. The user-provided parameter \emph{max_rearrangements} determines the two paths the function \emph{mbMCWtest} can follow to generate the collection of expected-by-chance mbMCW_BIs:
\itemize{
\item \emph{mbMCW exact testing}: If the number of distinct matched-paired measure rearrangements that can alter their initial set distribution is less than \emph{max_rearrangements}, the function \emph{mbMCWtest} calculates mbMCW_BIs for all possible data rearrangements.
\item \emph{mbMCW approximated testing}: If the number of distinct matched-paired measure rearrangements that can alter their initial set distribution is greater than \emph{max_rearrangements}, the function \emph{mbMCWtest} performs N = \emph{max_rearrangements} random measure rearrangements to calculate the collection of expected-by-chance mbMCW_BIs.
}
\item It calculates \emph{P~upper~} and \emph{P~lower~} values as the fraction of expected-by-chance mbMCW_BIs that are higher or equal to and lower or equal to the observed mbMCW_BIs, respectively.
}
}
}
\examples{

test_temp <- tempdir()
extdata_v <- system.file("extdata", "example_vertical_mbMCWtest_data.csv", package = "MCWtests")
file.copy(extdata_v, test_temp)
extdata_h <- system.file("extdata", "example_horizontal_mbMCWtest_data.csv", package = "MCWtests")
file.copy(extdata_h, test_temp)
# running mbMCWtest with an ideal vertical entry dataset
path_v <- file.path(test_temp, "example_vertical_mbMCWtest_data.csv")
mbMCWtest_vertical_results <- mbMCWtest(path_v, 10)
print(mbMCWtest_vertical_results)
# running mbMCWtest with an ideal horizontal entry dataset
path_h <- file.path(test_temp, "example_horizontal_mbMCWtest_data.csv")
mbMCWtest_horizontal_results <- mbMCWtest(path_h, 10)
print(mbMCWtest_horizontal_results)

rm(test_temp)

}


## Installation

You can install the development version of MCWtests from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# install MCWtests from github repository
devtools::install_github("diazcastillo/MCWtests")
# install MCWtests with vignettes from github repository
devtools::install_github("diazcastillo/MCWtests", build_vignettes = TRUE)
```

# Introduction to MCW testing

Originally, Monte Carlo-Wilcoxon (MCW) tests were designed to determined
whether the differences between two sets of data were significantly
biased in the same direction when compared with what it would be
expected by chance. MCW tests proceed by calculating sum-of-ranks-based
bias indexes, hence the reference to Frank Wilcoxon who invented the
non-parametric rank-sum and signed-rank tests, before and after
rearranging the dataset multiple times, hence the Monte Carlo reference
often associated to analytical strategies based on repeated random
sampling.

The *MCWtests* package encompasses the original MCW test and three
variations that differ in the data structures and the specific questions
they interrogate.

- The **matched-measures univariate MCW (muMCW)** **test**, the original
  MCW test, assesses whether one set of inherently matched-paired
  measures is significantly biased in the same direction. For instance,
  muMCW tests can be used to analyze bodyweights or transcript
  abundances determined at two different timepoints for the same set of
  mice.

- The **unmatched-measures MCW (uMCW) test** assesses whether two sets
  of unmatched measures and their heterogeneity are significantly biased
  in the same direction. For instance, uMCW tests can be used to analyze
  bodyweights or transcript abundances determined for two sets of mice
  that have been maintained in different conditions.

- The **matched-measures bivariate MCW (mbMCW) test** assesses whether
  two sets of inherently matched-paired measures are significantly
  differentially biased in the same direction. For instance, mbMCW tests
  can be used to analyze bodyweights or transcript abundances determined
  at two different timepoints for two sets of mice that have been
  exposed to different conditions.

- The **bias-measures MCW (bMCW) test** assesses whether a set of
  measures of bias for a quantitative trait between two conditions or a
  subset of these bias measures are themselves significantly biased in
  the same direction. For instance, bMCW tests can be used to analyze
  bias indexes obtained using other MCW tests or fold change for
  transcript abundances spanning the entire transcriptome or only for
  genes located in specific genomic regions from two sets of mice
  exposed to different conditions.

# MCW testing process

Although each MCW test examines distinct data structures to address
slightly different questions, all MCW tests share two fundamental steps:

1.  **To quantitatively determine the extent and direction of the bias
    of the measure under analysis**, MCW tests calculate a bias index
    (BI) by summing ranks and dividing these sums by the maximum
    possible value of the sums. Consequently, BIs range from 1 to -1
    when the measure under analysis is completely biased in each
    possible direction.

2.  **To determine the significance of the BIs calculated for the
    user-provided dataset (observed BIs)**, a collection of
    expected-by-chance BIs is generated by rearranging the original
    dataset multiple times and calculating BIs for each iteration.
    *P<sub>upper</sub>* and *P<sub>lower</sub>* values are calculated as
    the fractions of expected-by-chance BIs that have values higher or
    equal to and lower or equal to the observed BIs, respectively.

Each MCW test employs the user-provided parameter *max_rearrangements*
to follow two alternative paths.

- **MCW exact tests.** If the number of distinct rearrangements that can
  be generated from the dataset under analysis is less than
  *max_rearrangements*, MCW tests will actually generate all possible
  data rearrangements to create the collection of expected-by-chance
  BIs. In this case, *P<sub>upper</sub>* and *P<sub>lower</sub>* values
  will be exact estimations of the likelihood of obtaining BIs with
  equal or more extreme values compared to observed BIs with datasets of
  the same size and range but different internal structures.

- **MCW approximated tests.** If the number of distinct rearrangements
  that can be generated from the dataset under analysis is greater than
  *max_rearrangements*, MCW tests will perform a specified number of
  random data rearrangements, equal to the value of
  *max_rearrangements*, to generate the collection of expected-by-chance
  BIs. In this case, *P<sub>upper</sub>* and *P<sub>lower</sub>* values
  will represent approximate estimations of the likelihood of obtaining
  BIs with equal or more extreme values compared to observed BIs with
  datasets of the same size and range but different internal structures.

# Further reading

Please refer to the *MCWtests* vignette for comprehensive descriptions
of the package’s structure, the four MCW tests that can be executed, and
examples.

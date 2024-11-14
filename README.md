
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MR LOVA (Latent Outcome Variable Approach)

<!-- badges: start -->
<!-- badges: end -->

**`MRLOVA`** is an R package designed to enhance Mendelian Randomization
(MR) analysis by incorporating a latent phenotype of the outcome, which
enables a more precise disentanglement of horizontal and vertical
pleiotropy effects. This feature allows for a more explicit assessment
of the exclusion restriction assumption and iteratively refines causal
estimates using the expectation-maximization algorithm. Compared to
existing MR methods, this approach offers a unique and potentially more
accurate framework.

In addition to its core functionality, the package includes permutation
testing, providing a non-parametric method for hypothesis testing that
controls the type 1 error rate under the null hypothesis without relying
on specific distributional assumptions.

**`MRLOVA`** also features two additional statistical tests:

- **Directional Pleiotropy Test:** Assesses whether genetic variants
  influence the outcome through pathways other than the exposure of
  interest.
- **InSIDE Assumption Test:** Evaluates whether the assumption that
  Instrument Strength Independent of Direct Effect (InSIDE) holds.

There are three functions:

- **`mr_lova()`** main function that performs MR analysis and provides
  causal effect estimate.
- **`InSIDE()`** performs Instrument Strength Independent of Direct
  Effect (InSIDE) assumption test.
- **`directional()`** performs directional pleiotropy test.

To view the help pages for functions in this package, prepend the
function name with a question mark.

``` r
#library(MRLOVA)
#?mr_lova
#?InSIDE
#?directional
```

## Installation

You can install the development version of **MRLOVA** from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lamessad/MRLOVA")
```

## Steps for Mendelian Randomization Analysis with MRLOVA

**1. Load Summary Statistics:**

- Begin by loading the summary statistics for the instrumental
  variables.

**2. Harmonize Data:**

- Check the alleles and effects between the exposure and outcome
  studies.
- Flip the sign of `betaY` if the effect allele in the exposure study
  differs from the effect allele in the outcome study.

**3. Standardize data:**

- Standardize the summary statistics before performing the MR analysis.
  The summary statistics are standardized as z-statistics and then
  rescaled by the sample size for both continuous and binary traits.

**4. Perform MR Analysis:**

- Use the standardized and harmonized summary statistics to perform the
  Mendelian Randomization analysis.

**Note**: If you are using the **“TwoSampleMR”** package, you can
leverage the harmonized summary data to create the input parameters for
`MRLOVA` as follows:

- **Continuous outcome** - the summary statistics are standardized as
  z-statistics and then rescaled by the sample size.

- **Binary outcome** - If the effect size of the outcome GWAS summary
  data was reported in Odds Ratio (OR), convert the OR to log(OR). Then,
  the summary statistics are standardized as z-statistics and then
  rescaled by the sample size.

``` r
#library(TwoSampleMR)
#library(MRLOVA)
#dat=harmonise_data(exposure_data, outcome_data)
#detaX=dat$z.exposure/sqrt(dat$samplesize.exposure)
#betaY=dat$z.outcome/sqrt(dat$samplesize.outcome) 
#betaXse=1/sqrt(dat$samplesize.exposure)
#betaYse=1/sqrt(dat$samplesize.outcome)
#ny=dat$samplesize.outcome
#est = mr_lova(betaY, betaX, betaYse, betaXse, ny, permutn = 1000,log_file = "log.txt")  
```

## Example

Here is an example demonstrating how to apply MRLOVA methods to infer
the causal effect from exposure to outcome, assuming the first two steps
mentioned above have been carried out.

``` r
#MR analysis
library(MRLOVA)
data("dat")#not 
head(dat)
#>        betaX    betaXse        betaY    betaYse    ny
#> 1 0.85315066 0.01493736  1.103462683 0.02275630 1e+05
#> 2 0.08639250 0.01081406  0.070145019 0.01635072 1e+05
#> 3 0.06845758 0.01171665  0.058050158 0.01789296 1e+05
#> 4 0.56542496 0.01039496  0.294865158 0.01602819 1e+05
#> 5 0.03413034 0.01511759 -0.006703383 0.02312029 1e+05
#> 6 0.06143616 0.01089646  0.016728266 0.01663377 1e+05
dat$nx=dat$ny # sample size of exposure
betaX = dat$betaX/dat$betaXse/sqrt(dat$nx)
betaY = dat$betaY/dat$betaYse/sqrt(dat$ny)
betaXse = 1/sqrt(dat$nx)
betaYse = 1/sqrt(dat$ny)
ny = dat$ny
est = mr_lova(betaY, betaX, betaYse, betaXse, ny, permutn = 1000) 
#> Warning in mr_lova(betaY, betaX, betaYse, betaXse, ny, permutn = 1000): To get
#> a more precise p-value, it is recommended to increase the number of
#> permutations, given the p-value of causal effects = 1.17462227053672e-25
est
#> $CausEst
#> [1] 0.338072
#> 
#> $CausEstSE
#> [1] 0.01383628
#> 
#> $CausEstP
#> [1] 1.174622e-25
#> 
#> $IVs
#>   [1]  0.000000e+00  2.726007e-02  7.790426e-02  9.916483e-01  1.434340e-01
#>   [6]  2.109689e-01  8.345624e-68  1.281148e-04  5.792920e-01  4.060073e-01
#>  [11]  3.169787e-01  0.000000e+00  8.007780e-01  9.257839e-01  1.266553e-39
#>  [16]  8.154406e-01  4.525763e-38  6.405958e-04  7.576975e-01 4.858962e-234
#>  [21]  9.017435e-01  7.279847e-04  4.751810e-01  3.282271e-03  5.379763e-01
#>  [26]  1.700118e-01  9.558248e-01  6.680886e-01  1.016138e-02  4.648153e-03
#>  [31]  2.537098e-01  1.888917e-86  9.642804e-01  1.021426e-02  9.383683e-01
#>  [36]  5.601607e-18  2.234118e-01 5.278566e-311  9.622360e-01  8.652901e-01
#>  [41]  6.723178e-02  2.395406e-02  9.946421e-01  1.615649e-01  8.842703e-01
#>  [46]  0.000000e+00  4.322573e-01  1.252075e-01  9.039271e-01  9.074447e-01
#>  [51]  0.000000e+00  5.873736e-01  2.972120e-01  5.963484e-01  0.000000e+00
#>  [56]  6.681700e-01  0.000000e+00 1.141793e-140  3.080328e-01  2.109287e-01
#>  [61]  0.000000e+00  3.557337e-01 9.177200e-159  1.237713e-90  0.000000e+00
#>  [66]  3.565802e-01  2.648228e-02  6.289434e-01  1.262372e-01  1.673632e-03
#>  [71]  9.548446e-01  8.465089e-01  3.137614e-01  8.260025e-01 1.417906e-153
#>  [76]  1.459041e-02  4.677599e-01  1.133160e-01  0.000000e+00  1.174469e-11
#>  [81]  0.000000e+00  6.461811e-01  8.840330e-02  2.193040e-94  8.092369e-03
#>  [86]  6.103384e-03  0.000000e+00  0.000000e+00  2.546705e-01  1.408883e-01
#>  [91]  7.517597e-01  5.864770e-01  9.289636e-01  0.000000e+00  0.000000e+00
#>  [96]  1.759770e-01  9.665809e-01  0.000000e+00  6.080421e-01  4.277376e-02
#> 
#> $Valid
#>  [1]  3  4  6  9 11 14 16 19 25 26 27 33 35 37 39 40 41 43 45 47 48 49 50 52 53
#> [26] 54 59 60 69 72 73 74 77 78 83 89 91 92 96 97 99
#> 
#> $sig_v
#>           5% 
#> 1.694689e-13 
#> 
#> $corrected_p
#> [1] 0
```

#### Warnings

The package includes built-in warnings to help ensure the robustness of
your analysis. Please be aware of the following conditions and
corresponding warnings:

**Convergence Warnings** - If the number of iterations is less than 3 or
greater than 10, the package will issue a warning:

``` r
#warning("Please check convergence.")
#log_message("Warning: Please check convergence.", log_file)
```

This indicates potential issues with the convergence process. Details of
the convergence are saved in the `log.txt` file in the directory where
you are running the analysis. Reviewing this log file can provide
insights and help troubleshoot potential problems.

**Permutation Warnings** - Performing 1000 permutations is typically
sufficient for most analyses. However, for a more accurate p-value,
increasing the number of permutations is recommended. Ideally, the
number of permutations should exceed `1 / causal effect p-value`,
although this can be impractical. If the number of permutations is less
than `1/est$CausEstP` (where `est$CausEstP` is the estimated p-value
without permutation), the package will issue a warning:

``` r
#warning("#To get a more precise p-value, it is recommended to increase the number of permutations, given the p-value of causal effects  =  ", est$CausEstP)
#log_message(paste("# To get a more precise p-value, increasing the number of permutations is recommended, given the p-value of causal effects, p-value of causal effect  =  ", est$CausEstP), log_file)
```

Here are two outputs of the permutation analysis: `est$sig_v`(the 5th
percentile of the permuted p-values) and `corrected_p`(the corrected
p-value based on the permutation distribution). The permuted p-values
form an empirical null distribution under the assumption of no effect.
The 5th percentile of this distribution provides a threshold or critical
value for significance. If `est$CausEstP` is less than `est$sig_v`, it
indicates that the observed p-value is significantly smaller than what
would be expected under the null hypothesis. The `corrected_p` is the
proportion of permuted p-values that are less than or equal to the
observed p-value (`est$CausEstP`), representing the corrected p-value
based on permutation.

### InSIDE

**1. Using all instrumental variables:**

``` r
# InSIDE assumption test
inside_all = InSIDE(betaY, betaX, as.numeric(est$CausEst))  
inside_all
#> [[1]]
#>       cor 
#> 0.7927623 
#> 
#> $p_value
#> [1] 8.572826e-23
```

**2. Using only valid instrumental variables:** The `MRLOVA` package
identifies valid instrumental variables, and their indices are listed in
the analysis output as `est$Valid`.

``` r
# InSIDE assumption test
inside = InSIDE(betaY[est$Valid], betaX[est$Valid], as.numeric(est$CausEst))  
inside
#> [[1]]
#>         cor 
#> -0.01718479 
#> 
#> $p_value
#> [1] 0.9150736
```

### directional

**1. Using all instrumental variables:**

``` r
#directional pleiotropy test
dir_pleiotropy=directional(betaY, betaX, as.numeric(est$CausEst))
dir_pleiotropy
#> $estimate
#> [1] 0.02520747
#> 
#> $p_value
#> [1] 2.094721e-05
```

**2. Using only valid instrumental variables:**

``` r
#directional pleiotropy test
dir_pleiotropy=directional(betaY[est$Valid], betaX[est$Valid], as.numeric(est$CausEst))
dir_pleiotropy
#> $estimate
#> [1] 3.091524e-05
#> 
#> $p_value
#> [1] 0.9250098
```

**Note**: This novel statistical test has more statistical power than
the MR-Egger intercept for detecting directional pleiotropy, while also
controlling the Type I error rate(refer the manuscript).

## Contact

Please contact Lamessa Amente (<lamessa.amente@mymail.unisa.edu.au>) or
Hong Lee (<hong.lee@unisa.edu.au>) if you have any queries.

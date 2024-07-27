
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRLOVA

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

**`MRLOVA`** also features two additional tests:

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

**2. Harmonize Data:** - Check the allele frequencies of palindromic
SNPs. - Flip the sign of `betaY` if the effect allele in the exposure
study differs from the effect allele in the outcome study.

**3. Harmonize Data:** - Standardize the summary statistics before
performing the MR analysis. This approach follows the guidelines
suggested in the MRmix software (Qi and Chatterjee, Nat Commun 2019).

**4. Perform MR Analysis:** - Use the standardized and harmonized
summary statistics to perform the Mendelian Randomization analysis.

**Note**: If you are using the **“TwoSampleMR”** package, you can
leverage the harmonized summary data to create the input parameters for
`MRLOVA` as follows:

``` r
#library(TwoSampleMR)
#library(MRLOVA)
#dat=harmonize_data(exposure_data, outcome_data)
#detaX=dat$z.exposure/sqrt(dat$samplesize.exposure)
#betaY=dat$z.outcome/sqrt(dat$samplesize.outcome)
#betaXse=1/sqrt(dat$samplesize.exposure)
#betaYse=1/sqrt(dat$samplesize.outcome)
#ny=dat$samplesize.outcome
#est = mr_lova(betaY, betaX, betaYse, betaXse, ny, permutn = 1000,log_file = "log.txt")  
```

## Example

Here is an example demonstrating how to apply MRLOVA methods to infer
the causal effect from exposure to outcome, assuming the first three
steps mentioned above have been carried out. \### mr_lova

``` r
#MR analysis
library(MRLOVA)
data("dat")
head(dat)
#>        betaX    betaXse        betaY    betaYse    ny
#> 1 0.85315066 0.01493736  1.103462683 0.02275630 1e+05
#> 2 0.08639250 0.01081406  0.070145019 0.01635072 1e+05
#> 3 0.06845758 0.01171665  0.058050158 0.01789296 1e+05
#> 4 0.56542496 0.01039496  0.294865158 0.01602819 1e+05
#> 5 0.03413034 0.01511759 -0.006703383 0.02312029 1e+05
#> 6 0.06143616 0.01089646  0.016728266 0.01663377 1e+05
```

``` r
betaX = dat$betaX 
betaY = dat$betaY
betaXse = dat$betaXse
betaYse = dat$betaYse
ny = dat$ny
est = mr_lova(betaY, betaX, betaYse, betaXse, ny, permutn = 1000) 
#> Warning in mr_lova(betaY, betaX, betaYse, betaXse, ny, permutn = 1000): #
#> permutations may not be sufficient, given the p-value of causal effects, causal
#> p-value = 4.64046859466231e-28
```

``` r
est
#> $CausEst
#>  Estimate 
#> 0.5142296 
#> 
#> $CausEstSE
#>         SE 
#> 0.02053385 
#> 
#> $CausEstP
#>            P 
#> 4.640469e-28 
#> 
#> $IVs
#>   [1] 6.615865e-190  1.113516e-01  1.971081e-01  7.952652e-01  2.912525e-01
#>   [6]  3.659330e-01  7.348671e-37  5.409938e-03  6.850492e-01  5.468950e-01
#>  [11]  4.771970e-01  0.000000e+00  8.501856e-01  9.524966e-01  5.919295e-22
#>  [16]  8.590236e-01  3.261935e-21  1.260326e-02  8.201112e-01 1.028667e-126
#>  [21]  9.282667e-01  1.473534e-02  6.035096e-01  3.403586e-02  6.538054e-01
#>  [26]  3.183142e-01  9.628735e-01  7.583114e-01  6.196335e-02  4.003864e-02
#>  [31]  4.053961e-01  2.671107e-46  9.729690e-01  6.255681e-02  9.624724e-01
#>  [36]  4.075133e-10  3.791737e-01 8.348075e-166  9.744868e-01  8.942730e-01
#>  [41]  1.891448e-01  1.007077e-01  9.928422e-01  3.104539e-01  9.359081e-01
#>  [46] 3.226988e-265  5.661112e-01  2.595720e-01  9.283617e-01  9.240324e-01
#>  [51]  0.000000e+00  6.937718e-01  4.490627e-01  7.017550e-01  0.000000e+00
#>  [56]  7.562262e-01 3.934105e-286  9.721216e-76  4.577651e-01  3.586034e-01
#>  [61] 1.780989e-210  5.034448e-01  1.161396e-84  1.360985e-48  0.000000e+00
#>  [66]  5.047572e-01  1.039080e-01  7.184237e-01  2.726579e-01  2.271982e-02
#>  [71]  9.640004e-01  8.906684e-01  4.637285e-01  8.699707e-01  1.500034e-82
#>  [76]  7.399201e-02  5.970349e-01  2.482544e-01 1.556953e-215  8.173761e-07
#>  [81] 1.192874e-264  7.392510e-01  2.151794e-01  6.386109e-52  5.648529e-02
#>  [86]  4.667167e-02  0.000000e+00 6.266908e-210  4.030541e-01  2.866828e-01
#>  [91]  8.123083e-01  6.865125e-01  9.490333e-01  0.000000e+00  0.000000e+00
#>  [96]  3.254062e-01  9.655269e-01  0.000000e+00  7.134955e-01  1.408646e-01
#> 
#> $Valid
#>  [1]  2  3  4  6  9 11 14 16 19 25 26 27 33 34 35 37 39 40 41 42 43 45 47 48 49
#> [26] 50 52 53 54 59 60 67 69 72 73 74 77 78 83 85 89 91 92 96 97 99
#> 
#> $sig_v
#>           5% 
#> 1.314388e-13 
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

**Permutation Warnings** - If the number of permutations is less than
`1/est$CausEstP` of the estimate without permutation, the package will
issue a warning:

``` r
#warning("# permutations may not be sufficient, given the p-value of causal effects, causal p-value = ", mrlova_result$CausEstP)
#log_message(paste("# permutations may not be sufficient, given the p-value of causal effects, causal p-value =", mrlova_result$CausEstP), log_file)
```

### InSIDE

``` r
# InSIDE assumption test
inside = InSIDE(betaY, betaX, as.numeric(est$CausEst))  
inside
#> [[1]]
#>       cor 
#> 0.7847893 
#> 
#> $p_value
#> [1] 4.42478e-22
```

### directional

``` r
#directional pleiotropy test
dir_pleiotropy=directional(betaY, betaX, as.numeric(est$CausEst))
dir_pleiotropy
#> $estimate
#> [1] 0.1338896
#> 
#> $p_value
#> [1] 1.620467e-05
```

#' Mendelian Randomization Analysis
#'
#' @description This function conducts Mendelian randomization analysis based on the latent phenotype of the outcome by explicitly excluding vertical pleiotropy variants. It iteratively refines the causal relationship through an Expectation-Maximisation (EM) algorithm. MRLOVA takes GWAS summary statistics as inputs to estimate the causal effects of one trait on another. We recommend using summary statistics on a standardized scale.
#' @param betaY beta of the outcome, recommended to be in standardized scale.
#' @param betaX beta of the exposure, recommended to be in standardized scale.
#' @param betaYse standard error of the outcome, recommended to be in standardized scale.
#' @param betaXse standard error of the exposure, recommended to be in standardized scale.
#' @param ny sample size of outcome GWAS.
#' @param gwas_p p-value of association between variants and the latent phenotype ( the direct genetic effect), Default: 5e-2.
#' @param gwas_p2 p-value of exposure GWAS, Default: 5e-8.
#' @param permutn number of permutations, Default: 0.
#' @param log_file save the iteration to check the convergence, Default: NULL.
#' @details None.
#' @keywords Mendelian randomization.
#' @importFrom stats cor ecdf lm pchisq pt quantile
#' @export
#' @return A list that contains
#' \item{CausEst}{Estimate of causal effect.}
#' \item{CausEstSE}{Standard error of causal effect estimate.}
#' \item{CausEstP}{P-value from the z test for the causal effect.}
#' \item{IVs}{P-values for the direct causal effect of all Instrumental Variants (IV) on the outcome.}
#' \item{Valid}{Index of valid  instrumental variants.}
#' \item{corrected_p-value}{P-value from permutation tests.}



mr_lova <- function(betaY, betaX, betaYse, betaXse, ny, gwas_p = 5e-2, gwas_p2 = 5e-8, permutn = 0, log_file = NULL) {

  log_message <- function(message, log_file) {
    if (!is.null(log_file)) {
      cat(message, file = log_file, append = TRUE, sep = "\n")
    }
  }

  mrlova <- function(betaY, betaX, betaYse, betaXse, ny, gwas_p, gwas_p2, log_file) {
    sv1 <- seq_along(betaY)
    svv <- integer(0)
    ivw <- IVW(betaY[svv], betaX[svv], betaYse[svv])

    yi <- 0
    tau_t <- -99999
    a3_p <- pchisq((betaX / betaXse)^2, df = 1, lower.tail = FALSE)
    while (tau_t != ivw[1] & yi <= 10) {
      yi <- yi + 1
      tau_t <- ivw[1]
      gwas5 <- matrix(0, nrow = length(betaY), ncol = 4)
      pcor <- cor(betaY, betaX)
      gwas5[, 1] <- betaY - tau_t * betaX
      gwas5[, 2] <- sqrt(betaYse^2 + (tau_t^2 / ny) - (2 * tau_t * pcor / ny))
      gwas5[, 3] <- gwas5[, 1] / gwas5[, 2]
      gwas5[, 4] <- pchisq(gwas5[, 3]^2, df = 1, lower.tail = FALSE)
      svv <- sv1[gwas5[, 4] > gwas_p & a3_p < gwas_p2]

      log_message(paste("Iteration:", yi, "tau_t:", tau_t, "ivw[1]:", ivw[1]), log_file)

      ivw <- IVW(betaY[svv], betaX[svv], betaYse[svv])
    }

    if (yi < 3 | yi > 10) {
      warning("Please check convergence.")
      log_message("Warning: Please check convergence.", log_file)
    }

    ivw <- IVW(betaY[svv], betaX[svv], betaYse[svv])
    result <- list(CausEst = ivw[1], CausEstSE = ivw[2], CausEstP = ivw[4], SNPP = gwas5[, 4], Valid = svv)
    return(result)
  }

  mrlova_result <- mrlova(betaY, betaX, betaYse, betaXse, ny, gwas_p, gwas_p2, log_file)

  if (permutn > 0) {

    permutp <- numeric(permutn)

    for (zi in seq_len(permutn)) {
      sv2 <- sample(seq_along(betaX))
      mrlova_permutation <- mrlova(betaY[sv2], betaX, betaYse[sv2], betaXse, ny, gwas_p, gwas_p2, log_file)
      permutp[zi] <- mrlova_permutation$CausEstP
    }

    permutt <- quantile(permutp, probs = 0.05, na.rm = TRUE)
    permutt2 <- ecdf(permutp)(mrlova_result$CausEstP)

    if (1 / mrlova_result$CausEstP > permutn) {
      warning("# permutations may not be sufficient, given the p-value of causal effects, causal p-value = ", mrlova_result$CausEstP)
      log_message(paste("# permutations may not be sufficient, given the p-value of causal effects, causal p-value =", mrlova_result$CausEstP), log_file)
    }

    result <- list(CausEst=mrlova_result$CausEst,CausEstSE=mrlova_result$CausEstSE, CausEstP=mrlova_result$CausEstP, IVs=mrlova_result$SNPP, Valid=mrlova_result$Valid, sig_v=permutt, corrected_p=permutt2)
  } else {

    result <- mrlova_result
  }

  return(result)
}




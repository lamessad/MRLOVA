#' Instrument Strength Independent of Direct Effect (InSIDE) assumption
#'
#' @description This function calculates direct genetic effect of outcome  (`u`) from the difference between observed effect sizes of the outcome (`betaY`) and the predicted effect sizes based on the exposure (`betaX`) and a given causal estimate (`causEst`). It then performs a correlation test between these estimated the direct genetic effect of outcome  (`u`) and the exposure effect sizes.
#'
#' @param betaY Numeric vector of effect sizes for the outcome.
#' @param betaX Numeric vector of effect sizes for the exposure.
#' @param causEst Numeric value representing the estimated causal effect of the exposure on the outcome.
#'
#' @details The function calculates direct genetic effect of outcome `u` as `betaY - causEst * betaX`. Then performs a Pearson correlation test between the `u` and `betaX` to assess the relationship between the  exposure genetic effect and the direct genetic effect of the outcome. The function performs a pearson correlation test where the null hypothesis is that the correlation is zero.
#' @importFrom stats cor.test
#' @return A list containing:
#' \item{cor}{The estimate of the Pearson correlation coefficient between the direct genetic effect of the outcome and genetic effect of the exposure.}
#' \item{p_value}{The p-value from the correlation test.}
#'
#' @export
#'
InSIDE <- function(betaY, betaX, causEst) {
  u <- betaY - causEst * betaX
  cor_result <- cor.test(u, betaX)
  return(list(cor = cor_result$estimate, p_value = cor_result$p.value))
}

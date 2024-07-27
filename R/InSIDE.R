#' Instrument Strength Independent of Direct Effect (InSIDE) assumption test
#'
#' @description This function calculates direct genetic effect of outcome  (`u`) from the difference between observed effect sizes of the outcome (`betaY`) and the predicted effect sizes based on the exposure (`betaX`) and a given causal estimate (`causEst`). It then performs a correlation test between these estimated the direct genetic effect of outcome  (`u`) and the exposure effect sizes.
#'
#' @param betaY beta of the outcome, recommended to be in standardized scale.
#' @param betaX beta of the exposure, recommended to be in standardized scale.
#' @param causEst the estimated causal effect of the exposure on the outcome.
#'
#' @details The function calculates direct genetic effect of the outcome `u` as `betaY - causEst * betaX`. It then performs a Pearson correlation test between `u` and `betaX` to estimate the correlation between the  exposure genetic effect and the direct genetic effect of the outcome.
#' @importFrom stats cor.test
#' @return A list containing:
#' \item{cor}{The estimate of the Pearson correlation coefficient between the direct genetic effect of the outcome and the genetic effect of the exposure.}
#' \item{p_value}{The p-value from the Pearson correlation test.}
#'
#' @export
#'
InSIDE <- function(betaY, betaX, causEst) {
  u <- betaY - causEst * betaX
  cor_result <- cor.test(u, betaX)
  return(list(cor_result$estimate, p_value = cor_result$p.value))
}

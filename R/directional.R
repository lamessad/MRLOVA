#' Perform Directional Pleiotropy Test
#'
#' @description This function performs a one-sample t-test to determine whether the mean of the direct genetic effect on the outcome is significantly different from zero. This test helps evaluate the presence of directional pleiotropy by checking if the mean of the direct genetic effect on the outcome deviates significantly from zero.
#'
#' @param betaY beta of the outcome, recommended to be in standardized scale.
#' @param betaX beta of the exposure, recommended to be in standardized scale.
#' @param causEst the estimated causal effect of the exposure on the outcome.
#'
#' @details The function performs a t-test where the null hypothesis is that the mean of `u` is zero. The result includes the estimated mean of `u` and the p-value from the t-test.
#' @importFrom stats cor.test t.test
#' @return A list containing:
#' \item{mean}{The mean of the direct genetic effect of the outcome.}
#' \item{p_value}{The p-value from the t-test.}
#'
#' @export
directional <- function(betaY, betaX, causEst) {
  u <- betaY - causEst * betaX
  ttest_result <- t.test(u, mu = 0)
   return(list(estimate = as.numeric(ttest_result$estimate), p_value = ttest_result$p.value))
}

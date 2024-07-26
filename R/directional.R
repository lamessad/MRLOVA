#' Perform Directional Pleiotropy Test
#'
#' @description This function performs a one-sample t-test to assess whether the mean of the residuals (`u`) is significantly different from zero. This test helps in evaluating the presence of directional pleiotropy by checking if the residuals deviate significantly from zero.
#'
#' @param u Numeric vector of direct genetic effect of the outcome.
#'
#' @details The function performs a t-test where the null hypothesis is that the mean of `u` is zero. The result includes the estimated mean of `u` and the p-value from the t-test.
#' @importFrom stats cor.test t.test
#' @return A list containing:
#' \item{mean}{The mean of the direct genetic effect of outcome `u` estimated from the t-test.}
#' \item{p_value}{The p-value from the t-test.}
#'
#' @export
directional <- function(u) {
  # Perform t-test
  ttest_result <- t.test(u, mu = 0)

  # Return results
  return(list(mean = ttest_result$estimate, p_value = ttest_result$p.value))
}

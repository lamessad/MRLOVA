#' Inverse-Variance Weighted (IVW) Estimation
#'
#' @description This function performs Inverse-Variance Weighted (IVW) estimation of causal effect using linear regression. It estimates the causal effect of an exposure on an outcome by weighting the observations inversely proportional to their variance. The function returns the estimate of the causal effect, its standard error, t-value, and p-value.
#'
#' @param betaY Numeric vector of effect sizes for the outcome.
#' @param betaX Numeric vector of effect sizes for the exposure.
#' @param betaYse Numeric vector of standard errors for the outcome.
#'
#' @details The function uses linear regression to estimate the causal effect of the exposure on the outcome, with weights based on the inverse of the variance of the outcome's effect sizes. The standard error is adjusted based on the model's residual standard error. The p-value is computed from the t-statistic.
#'
#' @return A named vector with the following elements:
#' \item{Estimate}{The estimate of the causal effect.}
#' \item{SE}{The standard error of the estimate.}
#' \item{t value}{The t-statistic for the estimate.}
#' \item{P}{The p-value for the t-test.}
#'

IVW= function(betaY, betaX, betaYse) {

  beta=0;se=1;p=1
  if (length(betaY) > 0) {
    wv = 1/betaYse^2
    mod1 = lm(betaY ~ -1 + betaX, weights=wv)
    sod1 = summary(mod1)
    df   = length(betaY) - 1
    beta = sod1$coefficients[1,1]
    #se   = sod1$coefficients[1,2]/sod1$sigma
    se   = sod1$coefficients[1,2]/min(1,sod1$sigma) #check
    p    = 2*(pt(abs(beta/se),df, lower.tail=F))
  }
  z    = c(beta,se,beta/se,p)
  names(z)= c("Estimate","SE","t value","P")
  return(z)
}

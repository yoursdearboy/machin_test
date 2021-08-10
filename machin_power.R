#' Estimate power of test for comparing survival proportions at fixed time.
#' It's advised to use simulation, see README.
#'
#' @param s1 Survival proportion in the 1st group.
#' @param s2 Survival proportion in the 2nd group.
#' @param e2 Effective sample size at the fixed time in the 2nd group (num. still at risk / survival probability).
#' @param k Ratio between sample size e1 / e2.
#' @param mu A shift for non-inferiority / superiority testing.
#' @param sig.level Significance level \alpha.
#'
#' @return Statistical power of the test.
machin_power <- function(s1, s2, e2, k = 1, mu = 0, sig.level = .05, alternative = c("two.sided", "one.sided")) {
  alternative <- match.arg(alternative)

  e1 <- k * e2
  z <- (s1 - s2 - mu) / sqrt((s1 * (1 - s1) / e1) + (s2 * (1 - s2) / e2))

  if (alternative == "two.sided") {
    pnorm(z - qnorm(1 - sig.level / 2)) + pnorm(-z - qnorm(1 - sig.level / 2))
  } else {
    pnorm(abs(z) - qnorm(1 - sig.level))
  }
}

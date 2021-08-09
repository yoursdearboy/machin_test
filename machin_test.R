library(survival)

#' Survival data by month for 49 patients with Dukes's C colorectal cancer randomly
#' assigned to receive either y linolenic acid (1) or control (2) treatment.
#' 
#' @references
#' \cite{Machin, Gardner. Calculating confidence intervals for survival time analyses, 1988.}
ccancer <- data.frame(time = c(1,5,6,6,9,10,10,10,12,12,12,12,12,13,15,16,20,24,24,27,32,34,36,36,44,
                               3,6,6,6,6,8,8,12,12,12,15,16,18,18,20,22,24,28,28,28,30,30,33,42),
                      status = c(0,0,1,1,0,1,1,0,1,1,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,
                                 0,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,1,0,0,0,1,0,0,1),
                      group = c(rep(1, 25), rep(2, 24)))

#' Test for comparing survival proportions at fixed time.
#' 
#' @param t Fixed time to compare at.
#' @param alternative A character string specifying the alternative hypothesis.
#' @param mu A shift for non-inferiority / superiority testing.
#' @param conf.level Confidence level of the returned confidence interval.
#' @return List with statistics.
#' \itemize{
#'   \item d - Survival proportion difference.
#'   \item e - Effective sample size n' at time t.
#'   \item se - Standard error of survival difference.
#'   \item z - z-score.
#'   \item pval - p-value based on z-score.
#'   \item ci - confidence interval for difference in survival proportion.
#' }
#' 
#' @examples
#' fit <- survfit(Surv(time, status) ~ group, ccancer)
#' machin_test(fit, t = 24)
#' 
#' @references
#' \cite{Machin D, Gardner MJ. Calculating confidence intervals for survival time analyses. Br Med J 1988; 296:1369-1371}
machin_test <- function(fit, t, alternative = c("two.sided", "less", "greater"), mu = 0, conf.level = .95) {
  alternative <- match.arg(alternative)
  conf.alpha <- 1 - conf.level

  times <- sort(unique(c(fit$time, t)))
  res <- summary(fit, time = times)
  surv <- res$surv[times == t]
  n.risk <- res$n.risk[times == t]
  n.event <- res$n.event[times == t]

  d <- diff(surv)
  e <- (n.risk - n.event) / surv
  se <- sqrt(sum(surv * (1 - surv) / e))
  ci <- c(d + qnorm(conf.alpha / 2) * se, d + qnorm(1 - conf.alpha / 2) * se)
  z <- (d - mu) / se
  if (alternative == "two.sided") {
    pval <- (1 - pnorm(abs(z))) * 2
  } else {
    pval <- pnorm(z, lower.tail = (alternative == "less"))
  }
  list(d = d, e = e, se = se, z = z, pval = pval, ci = ci)
}

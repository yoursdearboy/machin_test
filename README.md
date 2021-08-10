# Test for comparing survival proportions at fixed time

Test propose by David Machin in [Machin D, Gardner MJ. Calculating confidence intervals for survival time analyses. Br Med J 1988; 296:1369-1371](https://www.bmj.com/content/296/6633/1369).

```r
library(survival)
source("https://raw.githubusercontent.com/yoursdearboy/machin_test/main/machin_test.R")

fit <- survfit(Surv(time, status) ~ group, ccancer)
plot(fit, mark.time = T, col = 2:3)
machin_test(fit, t = 24)
```

# Power of the test

## Example 1. Equality test

```r
library(survival)
source("https://raw.githubusercontent.com/yoursdearboy/machin_test/main/machin_power.R")

s1 <- 0.60 # Survival proportion in the 1st group.
s2 <- 0.76 # Survival proportion in the 2nd group.

# Effective sample size at the fixed time in the 2nd group (num. still at risk / survival probability).
e2 <- 117
# In case of known censoring distribution it may be estimated at the fixed time (e.g. 3 years) like this
# (example is uniform censoring with min 1 year follow-up and 3 year accrual).
n2 <- 350
e2 <- n2 * (1 - punif(3, 1, 4))

machin_power(s1 = s1, s2 = s2, e2 = e2, k = 1, sig.level = 0.05)
```

or using a simulation

```r
nsim <- 1e4 # number of simulations
mint <- 1 # min follow-up
acct <- 3 # accrual time
maxt <- mint + acct
lambda1 <- 0.09 # hazard rate in the 1st group
lambda2 <- 0.09 * 1.9 # hazard rate in the 2nd group (1.9 is hazard ratio)
n <- 350 # sample size
tt <- 3 # fixed time to test
alpha <- 0.05

sim <- replicate(nsim, {
  gr <- c(rep(1, n), rep(2, n))
  y <- c(rexp(n, lambda1), rexp(n, lambda2))
  cens <- runif(n * 2, mint, maxt)
  s <- as.numeric(y <= cens)
  t <- pmin(y, cens)
  fit <- survfit(Surv(t, s) ~ gr)
  machin_test(fit, tt)$p <= alpha
})

mean(sim) # power of the test
```

## Example 2. Non-inferiority test

```r
library(survival)
source("https://raw.githubusercontent.com/yoursdearboy/machin_test/main/machin_power.R")

s1 <- 0.72 # Survival proportion in the 1st group.
s2 <- 0.76 # Survival proportion in the 2nd group.
mu <- -0.1 # Non-inferiority margin.

# Effective sample size at the fixed time in the 2nd group (num. still at risk / survival probability).
e2 <- 117
# In case of known censoring distribution it may be estimated at the fixed time (e.g. 3 years) like this
# (example is uniform censoring with min 1 year follow-up and 3 year accrual).
n2 <- 350
e2 <- n2 * (1 - punif(3, 1, 4))

machin_power(s1 = s1, s2 = s2, e2 = e2, mu = -0.1, alternative = "one.sided", k = 1, sig.level = 0.05)
```

or using a simulation

```r
nsim <- 1e4 # number of simulations
mint <- 1 # min follow-up
acct <- 3 # accrual time
maxt <- mint + acct
lambda1 <- 0.09 # hazard rate in the 1st group
lambda2 <- 0.09 * 1.2 # hazard rate in the 2nd group (1.9 is hazard ratio)
mu <- -0.1 # non-inferiority margin
n <- 350 # sample size
tt <- 3 # fixed time to test
alpha <- 0.05

sim <- replicate(nsim, {
  gr <- c(rep(1, n), rep(2, n))
  y <- c(rexp(n, lambda1), rexp(n, lambda2))
  cens <- runif(n * 2, mint, maxt)
  s <- as.numeric(y <= cens)
  t <- pmin(y, cens)
  fit <- survfit(Surv(t, s) ~ gr)

  # using one-sided test
  machin_test(fit, tt, mu = mu, alternative = "greater")$p <= alpha

  # using confidence interval
  # machin_test(fit, tt, conf.level = 0.9)$ci[1] > mu
})

mean(sim) # power of the test
```

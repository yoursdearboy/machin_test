machin_power <- function(alpha = .05, s1, s2, e2, k = 1) {
  e1 <- k * e2
  z <- (s1 - s2) / sqrt((s1 * (1 - s1) / e1) + (s2 * (1 - s2) / e2))
  pnorm(z - qnorm(1 - alpha/2)) + pnorm(-z - qnorm(1 - alpha / 2))
}

machin_sample_size <- function(alpha = .05, beta = .2, t, tmin, tmax, s1, s2, k = 1) {
  d <- s1 - s2
  se <- (s1 * (1 - s1) / k) + (s2 * (1 - s2))

  e2 <- (((qnorm(1 - beta) + qnorm(1 - alpha / 2)) / d)^2) * se
  e1 <- k * e2

  r2 <- e2 * s2
  r1 <- e1 * s1

  n2 <- e2 / (1 - punif(t, tmin, tmax))
  n1 <- e2 / (1 - punif(t, tmin, tmax)) * k

  list(e1 = e1, e2 = e2, r1 = r1, r2 = r2, n1 = n1, n2 = n2)
}

# Demo
lambda <- .09
hr <- 1.9
t <- 3
tmin <- 1
tmax <- 4
n <- 700

machin_power(s1 = 1 - pexp(t, lambda),
             s2 = 1 - pexp(3, lambda * hr),
             e2 = (n / 2) * (1 - punif(t, tmin, tmax)),
             k = 1)

# or using simulation
rcontrol <- function(n) rexp(n, sim_lambda * sim_hr)
rinterv <- function(n) rexp(n, sim_lambda)
rcens <- function(n) runif(n, sim_tmin, sim_tmax)

Hmisc::spower(rcontrol = rcontrol,
              rinterv = rinterv,
              rcens = rcens,
              nc = n / 2, ni = n / 2,
              test = function(S, group) {
                fit <- survfit(Surv(S[,1], S[,2]) ~ group)
                z <- machin_test(fit, t)$z
                z^2
              },
              nsim = 1e3,
              alpha = .05)

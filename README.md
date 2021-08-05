# Test for comparing survival proportions at fixed time

Test propose by David Machin in [Machin D, Gardner MJ. Calculating confidence intervals for survival time analyses. Br Med J 1988; 296:1369-1371](https://www.bmj.com/content/296/6633/1369).

```r
library(survival)
source("https://raw.githubusercontent.com/yoursdearboy/machin_test/main/machin_test.R")

fit <- survfit(Surv(time, status) ~ group, ccancer)
plot(fit, mark.time = T, col = 2:3)
machin_test(fit, t = 24)
```

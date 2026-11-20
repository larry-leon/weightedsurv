#!/usr/bin/env Rscript
# calibrate.R -- TASK_ps_weight_results_by_time.md, section 2: propensity model logit e = a + b X for the cells
# Es9-400 and Es9-400-null. From the weights alone on one large draw of X, using the true propensity; no survival
# fitting. Population ESS ratio of arm a = (sum of the arm's weights)^2 / (n_a * sum of the arm's squared weights),
# with each subject entering arm a with probability e (treated) or 1 - e (control), i.e. in expectation:
#   treated: 1 / (P * E[1/e]),   control: 1 / ((1 - P) * E[1/(1 - e)]),   P = E[e].
# (The ratio does not depend on the weight type: stabilized and ipw weights differ by a constant within an arm.)
#
#   Rscript --vanilla dev/audit/ps_by_time/calibrate.R > dev/audit/ps_by_time/calibration.txt

set.seed(20261120)
X <- rnorm(1e6)
pop <- function(a, b) {
  e <- plogis(a + b * X); P <- mean(e)
  c(P = P, ess_treated = 1 / (P * mean(1 / e)), ess_control = 1 / ((1 - P) * mean(1 / (1 - e))))
}
ref <- pop(-0.5, 1)                                   # cells 31 / 32 (logit e = -0.5 + X)
cat(sprintf("draw: X ~ N(0, 1), n = 1e6, seed 20261120\n"))
cat(sprintf("cells 31/32, logit e = -0.5 + X: P(treated) = %.4f, population ESS ratio treated = %.4f, control = %.4f\n",
            ref["P"], ref["ess_treated"], ref["ess_control"]))

# a(b): P(treated) equal to that of cells 31/32
a_of_b <- function(b) uniroot(function(a) pop(a, b)["P"] - ref["P"], c(-5, 5), tol = 1e-10)$root
target <- 0.91
f <- function(b) pop(a_of_b(b), b)["ess_treated"] - target
b <- uniroot(f, c(0.05, 1), tol = 1e-8)$root
a <- a_of_b(b)
# rounded to two decimals for the run
b2 <- round(b, 2); a2 <- round(a_of_b(b2), 2)
v <- pop(a2, b2)
cat(sprintf("solution for treated ESS ratio %.2f: a = %.6f, b = %.6f\n", target, a, b))
cat(sprintf("fixed for the run (two decimals): a = %.2f, b = %.2f\n", a2, b2))
cat(sprintf("  P(treated) = %.4f (cells 31/32: %.4f; difference %.4f)\n", v["P"], ref["P"], v["P"] - ref["P"]))
cat(sprintf("  population ESS ratio treated = %.4f, control = %.4f\n", v["ess_treated"], v["ess_control"]))
cat(sprintf("PARAMS a=%.2f b=%.2f\n", a2, b2))

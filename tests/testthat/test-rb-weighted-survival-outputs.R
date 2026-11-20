# Rb: with subject weights whose positive values are not all equal, logrank_results is the robust score
# test at 0 of the weighted Cox model and cox_results$score its robust score test; cox_results$se added.

test_that("Rb: logrank_results and cox_results$score equal the coxph robust score tests (ties, strata)", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(300, 4, ties = ties); d$two <- 2
    dc  <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w"))
    dcs <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", strata.name = "s"))
    cpb <- survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow")
    cps <- survival::coxph(survival::Surv(time, status) ~ treat + survival::strata(s), data = d, weights = w, robust = TRUE, ties = "breslow")
    cpe <- survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE)
    expect_lt(abs(dc$logrank_results$chisq / cpb$rscore - 1), 1e-8)
    expect_lt(abs(dc$logrank_results$chisq / dc$z.score^2 - 1), 1e-8)
    expect_lt(abs(dcs$logrank_results$chisq / cps$rscore - 1), 1e-8)
    expect_lt(abs(dcs$logrank_results$chisq / dcs$z.score_stratified^2 - 1), 1e-8)
    expect_lt(abs(dc$cox_results$score / summary(cpe)$robscore[1] - 1), 1e-8)
    expect_lt(abs(dc$cox_results$se / sqrt(cpe$var[1, 1]) - 1), 1e-12)
    expect_identical(dc$logrank_results$method, "coxph_robust_score")
    expect_identical(dc$logrank_results$pval, 1 - stats::pchisq(dc$logrank_results$chisq, 1))
    expect_match(dc$logrank_results$logrank_text, "^Logrank p = ")
    # hr_ci and pval rest on cox_results$se
    expect_equal(unname(as.numeric(dc$cox_results$hr_ci)), unname(exp(log(dc$cox_results$hr) + c(-1, 1) * stats::qnorm(0.975) * dc$cox_results$se)), tolerance = 1e-10)
    # the scheme does not change logrank_results with subject weights (the rho = 0 test)
    d1 <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", scheme = "fh", scheme_params = list(rho = 1, gamma = 0)))
    expect_identical(d1$logrank_results$chisq, dc$logrank_results$chisq)
  }
})

test_that("Rb: check_results() agrees and does not warn on a weighted FH(0,0) result", {
  d <- ws_sim_ipw(300, 4)
  dc <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w"))
  expect_no_warning(r <- ws_capture(check_results(dc)))
  r <- check_results(dc, verbose = FALSE)
  expect_lt((max(r) - min(r)) / mean(c(max(r), min(r))), 1e-8)
})

test_that("Rb: no weights, equal weights and 0/1 weights keep survdiff and the model-based score test", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(300, 4, ties = ties); d$two <- 2; d$w01 <- rep(c(1, 1, 0), length.out = nrow(d))
    u  <- ws_quiet(df_counting(d, "time", "status", "treat"))
    c2 <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "two"))
    z  <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w01"))
    sd <- survival::survdiff(survival::Surv(time, status) ~ treat, data = d)
    expect_lt(abs(u$logrank_results$chisq / sd$chisq - 1), 1e-10)
    expect_identical(u$logrank_results$chisq, c2$logrank_results$chisq)
    expect_identical(u$logrank_results$method, "survdiff"); expect_identical(c2$logrank_results$method, "survdiff")
    expect_identical(z$logrank_results$method, "survdiff")
    cp <- survival::coxph(survival::Surv(time, status) ~ treat, data = d, robust = TRUE)
    expect_identical(unname(u$cox_results$score), unname(summary(cp)$sctest[1]))
    expect_lt(abs(u$cox_results$se / sqrt(cp$var[1, 1]) - 1), 1e-12)
  }
})

test_that("Rb: the NA placeholder carries se = NA", {
  d <- ws_sim(60, 3); d$status[d$treat == 1] <- 0; d$w <- stats::runif(nrow(d), 0.5, 2)
  r <- suppressWarnings(df_counting(d, "time", "status", "treat", weight.name = "w"))
  expect_true(is.na(r$cox_results$se))
  expect_true(is.finite(r$logrank_results$chisq))
})

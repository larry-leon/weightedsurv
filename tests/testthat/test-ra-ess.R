# Ra: df_counting()$ess, the effective sample size by arm.

test_that("Ra: ess equals the direct computation, in the order of quantile_results", {
  d <- ws_sim_ipw(300, 4); d$two <- 2
  ess <- function(w) sum(w)^2 / sum(w^2)
  dc <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", arms = c("Exp", "Ctl")))
  expect_identical(names(dc$ess), c("group", "n", "ess", "ratio"))
  expect_identical(dc$ess$group, c("Exp", "Ctl"))
  expect_identical(dc$ess$group, as.character(dc$quantile_results$group))
  expect_equal(dc$ess$n, c(sum(d$treat == 1), sum(d$treat == 0)))
  expect_lt(abs(dc$ess$ess[1] / ess(d$w[d$treat == 1]) - 1), 1e-10)
  expect_lt(abs(dc$ess$ess[2] / ess(d$w[d$treat == 0]) - 1), 1e-10)
  expect_equal(dc$ess$ratio, dc$ess$ess / dc$ess$n, tolerance = 1e-12)
  expect_true(all(dc$ess$ratio <= 1 + 1e-12))
})

test_that("Ra: ess = n without weights and with constant weights; zero-weight rows are not counted", {
  d <- ws_sim_ipw(300, 4); d$two <- 2; d$wz <- ifelse(seq_len(nrow(d)) %% 3 == 0, 0, d$w)
  u <- ws_quiet(df_counting(d, "time", "status", "treat")); c2 <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "two"))
  expect_lt(max(abs(u$ess$ess - u$ess$n)), 1e-9); expect_lt(max(abs(c2$ess$ess - c2$ess$n)), 1e-9)
  expect_equal(u$ess$ratio, c(1, 1), tolerance = 1e-12)
  z <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "wz"))
  p <- ws_quiet(df_counting(d[d$wz > 0, ], "time", "status", "treat", weight.name = "wz"))
  expect_equal(z$ess, p$ess)
  expect_equal(z$ess$n, c(sum(d$treat == 1 & d$wz > 0), sum(d$treat == 0 & d$wz > 0)))
})

test_that("Ra: nothing else in the result changes and no warning depends on ess", {
  d <- ws_sim_ipw(200, 9); d$wx <- ifelse(seq_len(nrow(d)) == 1, 1e4, 1)     # one subject carries most of the weight
  expect_no_warning(r <- ws_capture(df_counting(d, "time", "status", "treat", weight.name = "wx")))
  expect_lt(min(r$ess$ratio), 0.05)
  expect_true("ess" %in% names(r))
})

# B3: km_quantile() / kmq_calculations() / df_counting() quantiles agree with survfit.

sf_quant <- function(time, status, level, cl = 0.95, w = NULL) {
  fit <- if (is.null(w)) survival::survfit(survival::Surv(time, status) ~ 1, conf.int = cl, conf.type = "log")
         else survival::survfit(survival::Surv(time, status) ~ 1, weights = w, conf.int = cl, conf.type = "log")
  q <- stats::quantile(fit, probs = 1 - level)
  c(unname(q$quantile), unname(q$lower), unname(q$upper))
}
same_q <- function(a, b) all(is.na(a) == is.na(b)) && all(abs(a - b)[!is.na(a)] < 1e-10)

test_that("exact plateaus at the level: midpoint as in survfit (n = 2..60, five levels)", {
  for (lv in c(0.9, 0.75, 0.5, 0.25, 0.1)) for (n in seq(2, 60, by = 2)) {
    km <- KM_estimates(n:1, 1:n)
    expect_equal(km_quantile(1:n, km$S_KM, qprob = lv)$qhat, sf_quant(1:n, rep(1, n), lv)[1], tolerance = 1e-10)
  }
})

test_that("point estimate and unweighted limits equal survfit (ties, censoring, non-monotone limit curves)", {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  bad <- 0; n_cmp <- 0
  for (s in 1:60) {
    set.seed(5000 + s)
    n <- sample(c(10, 30, 60), 1)
    tt <- stats::rexp(2 * n, 0.08); cc <- stats::rexp(2 * n, 0.08 * c(0.01, 0.4, 2)[s %% 3 + 1])
    d <- data.frame(time = pmin(tt, cc), status = as.numeric(tt <= cc), treat = rep(0:1, each = n))
    if (s %% 2 == 0) d$time <- ceiling(d$time)
    k <- KM_diff(d, "time", "status", "treat")
    tf <- survival::aeqSurv(survival::Surv(d$time, d$status))[, "time"]
    i0 <- d$treat == 0
    last0 <- max(tf[i0]); sel <- k$at_points <= last0
    for (cl in c(0.8, 0.95)) for (lv in c(0.75, 0.5, 0.25)) {
      got <- km_quantile(k$at_points[sel], k$surv0[sel], sqrt(k$sig2_surv0[sel]), qprob = lv, conf_level = cl)
      n_cmp <- n_cmp + 1
      bad <- bad + !same_q(c(got$qhat, got$lower, got$upper), sf_quant(tf[i0], d$status[i0], lv, cl))
    }
  }
  expect_equal(bad, 0, info = sprintf("%d comparisons", n_cmp))
})

test_that("weighted point estimate equals survfit(weights =)", {
  d <- ws_sim(200, 41, ties = TRUE, weighted = TRUE)
  r <- suppressWarnings(df_counting(d, "time", "status", "treat", weight.name = "w"))
  for (g in 0:1) {
    dd <- d[d$treat == g, ]
    ref <- sf_quant(dd$time, dd$status, 0.5, w = dd$w)[1]
    expect_equal(r$quantile_results$quantile[if (g == 1) 1 else 2], ref, tolerance = 1e-10)
  }
})

test_that("plateau at the level running to the end of follow-up, and a curve that never reaches it", {
  d <- data.frame(time = c(1, 2, 3, 4, 5, 6, 1.5, 2.5, 3.5, 4.5), status = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0), treat = rep(0:1, c(6, 4)))
  r <- df_counting(d, "time", "status", "treat", get.cox = FALSE)
  # control: S = 4/6 from t = 2, never reaches 0.5 -> NA; treated: S = 0.75 from 1.5 to its end (4.5 < 6)
  expect_true(is.na(r$quantile_results$quantile[2]))
  r2 <- df_counting(d, "time", "status", "treat", qprob = 0.75, get.cox = FALSE)
  dd <- d[d$treat == 1, ]
  expect_equal(r2$quantile_results$quantile[1], sf_quant(dd$time, dd$status, 0.75)[1], tolerance = 1e-10)
  expect_equal(r2$quantile_results$quantile[1], (1.5 + 4.5) / 2)
})

test_that("df_counting(): quantile_check is survfit at probs = 1 - qprob and conf_level; no warning on valid data", {
  d <- ws_sim(200, 42)
  for (qp in c(0.25, 0.5, 0.75)) for (cl in c(0.9, 0.95)) {
    expect_no_warning(r <- df_counting(d, "time", "status", "treat", qprob = qp, conf_level = cl))
    dd <- d[d$treat == 0, ]
    expect_equal(unname(unlist(r$quantile_check[1, c("time", "lower", "upper")])), sf_quant(dd$time, dd$status, qp, cl), tolerance = 1e-10)
  }
})

test_that("df_counting(): a real quantile discrepancy is a warning whatever verbose is", {
  d <- ws_sim(120, 43)
  local_mocked_bindings(km_quantile = function(...) list(qhat = 999, lower = NA, upper = NA))
  expect_warning(df_counting(d, "time", "status", "treat", verbose = FALSE), "discrepancy in quantile")
  expect_error(df_counting(d, "time", "status", "treat", stop.onerror = TRUE), "discrepancy in quantile")
})

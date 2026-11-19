# B4: input validation and degenerate data.

null_device <- function(env = parent.frame()) {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  withr_defer <- function() if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev)
  do.call(on.exit, list(substitute(f(), list(f = withr_defer)), add = TRUE), envir = env)
}
n_warn <- function(expr) { k <- 0L; withCallingHandlers(expr, warning = function(w) { k <<- k + 1L; invokeRestart("muffleWarning") }); k }
kmd  <- function(d, ...) KM_diff(d, "time", "status", "treat", ...)
dfc  <- function(d, ...) df_counting(d, "time", "status", "treat", ...)
pkb  <- function(d, ...) { null_device(); plotKM.band_subgroups(d, "time", "status", "treat", draws.band = 0, ...) }

test_that("NA in time, event, treatment or weight is an error naming the column, in all three functions", {
  for (col in c("time", "status", "treat", "w")) {
    d <- ws_sim(80, 5, weighted = TRUE); d[[col]][4] <- NA
    expect_error(kmd(d, weight.name = "w"), sprintf("'%s'.*missing", col))
    expect_error(dfc(d, weight.name = "w"), sprintf("'%s'.*missing", col))
    expect_error(pkb(d, weight.name = "w"), sprintf("'%s'.*missing", col))
  }
  d <- ws_sim(80, 5); d$s <- rep(c("a", "b"), 40); d$s[3] <- NA
  expect_error(dfc(d, strata.name = "s"), "'s'.*missing")
})

test_that("non-finite or negative times are errors; zero is valid", {
  d <- ws_sim(80, 5)
  d1 <- d; d1$time[2] <- Inf
  d2 <- d; d2$time[2] <- -1
  d3 <- d; d3$time[2] <- 0
  for (f in list(kmd, dfc, pkb)) {
    expect_error(f(d1), "non-finite")
    expect_error(f(d2), "negative")
  }
  expect_silent(kmd(d3))
  expect_no_error(suppressWarnings(dfc(d3)))
})

test_that("negative or non-finite weights are errors", {
  d <- ws_sim(80, 5, weighted = TRUE); d$w[1] <- -1
  expect_error(kmd(d, weight.name = "w"), "'w'.*negative")
  expect_error(dfc(d, weight.name = "w"), "'w'.*negative")
  expect_error(pkb(d, weight.name = "w"), "'w'.*negative")
  d$w[1] <- Inf
  expect_error(kmd(d, weight.name = "w"), "'w'.*non-finite")
})

test_that("zero weights: KM as survfit on positive rows; df_counting self-checks silent, Cox and log-rank on positive rows", {
  d <- ws_sim(120, 6, weighted = TRUE); d$w[d$treat == 0][1:10] <- 0
  k <- kmd(d, weight.name = "w")
  dp <- d[d$w > 0, ]
  for (g in 0:1) {
    gg <- dp[dp$treat == g, ]
    ref <- summary(survival::survfit(survival::Surv(time, status) ~ 1, data = gg, weights = w), times = k$at_points, extend = TRUE)$surv
    expect_equal(if (g == 0) k$surv0 else k$surv1, ref, tolerance = 1e-10)
  }
  ws <- character(0)
  r <- withCallingHandlers(dfc(d, weight.name = "w"), warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("Discrepancy", ws)))
  expect_length(r$time, nrow(d))
  expect_length(r$w_hat, nrow(d))
  cx <- survival::coxph(survival::Surv(time, status) ~ treat, data = dp, weights = w, robust = TRUE)
  expect_equal(unname(r$cox_results$hr), unname(exp(stats::coef(cx))), tolerance = 1e-10)
  expect_equal(r$logrank_results$chisq, survival::survdiff(survival::Surv(time, status) ~ treat, data = dp)$chisq, tolerance = 1e-10)
})

test_that("event coded 0/1, logical and 1/2 give the same KM in all three functions; other codes are errors", {
  d <- ws_sim(100, 7)
  dl <- d; dl$status <- dl$status == 1
  d12 <- d; d12$status <- d12$status + 1
  ref_k <- kmd(d); ref_c <- suppressWarnings(dfc(d))
  ref_p <- pkb(d)
  for (dd in list(dl, d12)) {
    expect_equal(kmd(dd)$surv0, ref_k$surv0)
    expect_equal(kmd(dd)$surv1, ref_k$surv1)
    r <- suppressWarnings(dfc(dd))
    expect_equal(r$surv0, ref_c$surv0); expect_equal(r$z.score, ref_c$z.score)
    expect_equal(unname(r$cox_results$hr), unname(ref_c$cox_results$hr))
    expect_equal(pkb(dd)$mean, ref_p$mean)
  }
  d3 <- d; d3$status[1] <- 3
  for (f in list(kmd, dfc, pkb)) expect_error(f(d3), "'status'.*coded")
})

test_that("a missing column, including strata.name and weight.name, is an error naming it", {
  d <- ws_sim(60, 8)
  expect_error(dfc(d, strata.name = "nope"), "nope")
  expect_error(kmd(d, weight.name = "nope"), "nope")
  expect_error(pkb(d, weight.name = "nope"), "nope")
  expect_error(KM_diff(d, "time", "status", "arm"), "arm")
})

test_that("an empty arm is an error naming the arm (and the subgroup in plotKM.band_subgroups)", {
  d <- ws_sim(60, 9); d1 <- d; d1$treat <- 1
  expect_error(kmd(d1), "control arm")
  expect_error(dfc(d1), "control arm")
  expect_error(pkb(d1), "control arm")
  d$grp <- ifelse(d$treat == 1 & d$time > 5, "x", "y")
  d$grp[d$treat == 0] <- "y"
  expect_error(pkb(d, sg_labels = "grp == 'x'", sg_colors = "red"), "Subgroup 'grp == 'x''.*control arm")
})

test_that("stratified df_counting(): a stratum with one arm empty contributes nothing", {
  d <- ws_sim(150, 10); d$s <- rep(c("a", "b", "c"), 50)
  d$s[d$treat == 1 & d$s == "c"] <- "b"          # stratum c has no treated patient
  r <- suppressWarnings(dfc(d, strata.name = "s"))
  r2 <- suppressWarnings(dfc(d[d$s != "c", ], strata.name = "s"))
  expect_equal(r$lr_stratified, r2$lr_stratified, tolerance = 1e-10)
  expect_equal(r$sig2_lr_stratified, r2$sig2_lr_stratified, tolerance = 1e-10)
})

test_that("no events in an arm, all censored, n = 1 per arm: KM_diff returns, df_counting returns with one warning", {
  base <- ws_sim(60, 11)
  cases <- list(
    no_events_arm1 = transform(base, status = ifelse(treat == 1, 0, status)),
    all_censored = transform(base, status = 0),
    n1_per_arm = data.frame(time = c(3, 5), status = c(0, 1), treat = c(0, 1)),
    n1_tied = data.frame(time = c(8, 8), status = c(1, 1), treat = c(0, 1))
  )
  for (nm in names(cases)) {
    d <- cases[[nm]]
    expect_no_error(kmd(d))
    k <- n_warn(r <- dfc(d))
    expect_lte(k, 1)
    expect_true(is.list(r), info = nm)
    if (k == 1) expect_warning(dfc(d), "could not be computed")
  }
  r <- suppressWarnings(dfc(cases$no_events_arm1))
  expect_true(is.na(r$cox_results$hr))
  expect_true(all(r$surv1 == 1))
})

test_that("plotKM.band_subgroups(): no events in an arm without Maxtau is an error that mentions Maxtau", {
  d <- ws_sim(60, 12); d$status[d$treat == 0] <- 0
  expect_error(pkb(d), "max_tau is undefined.*Maxtau")
})

test_that("cox_rhogamma(): no events in an arm gives a warning, value unchanged", {
  d <- ws_sim(80, 13); d$status[d$treat == 1] <- 0
  dc <- suppressWarnings(dfc(d, get.cox = FALSE))
  expect_warning(f <- cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0)), "not identified")
  wt <- rep(1, length(dc$time))
  root <- find_cox_root(dc$time, dc$delta, dc$z, dc$w_hat, wt)$root
  expect_equal(f$fit$bhat, root)
})

test_that("KM_diff(draws.band < 20, show_resamples = TRUE) draws what there is", {
  null_device()
  for (b in c(1, 5, 19)) expect_no_error(kmd(ws_sim(80, 14), draws.band = b, show_resamples = TRUE))
})

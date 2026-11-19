#!/usr/bin/env Rscript
# static_leads_check.R -- small executable confirmations for static-read leads S1-S9
# and for new findings N2, N3, N5 (origin lines).  Local tree via pkgload::load_all.
# Run from the repo root:  Rscript dev/audit/static_leads_check.R > dev/audit/static_leads_check.txt 2>&1
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
grDevices::pdf(NULL)
run <- function(expr) {
  ws <- character(0)
  tc <- textConnection("out", "w", local = TRUE); sink(tc); on.exit({ sink(); close(tc) }, add = TRUE)
  val <- tryCatch(withCallingHandlers(expr,
           warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }),
         error = function(e) structure(conditionMessage(e), class = "err"))
  sink(); close(tc); on.exit()
  list(value = if (inherits(val, "err")) NULL else val, error = if (inherits(val, "err")) unclass(val) else NULL, warnings = ws, stdout = out)
}
where <- function(expr) {   # innermost package frame of an error
  r <- tryCatch(withCallingHandlers(expr, error = function(e) {
    calls <- sys.calls(); txt <- vapply(calls, function(x) paste(deparse(x, nlines = 1), collapse = ""), "")
    assign(".last_calls", utils::tail(txt[-length(txt)], 6), envir = .GlobalEnv)
  }), error = function(e) conditionMessage(e))
  r
}
say <- function(id, ...) cat(sprintf("[%s] ", id), ..., "\n", sep = "")
set.seed(1)
n <- 120; d <- data.frame(treat = rep(0:1, each = n / 2)); d$time <- rexp(n, 0.08); cc <- runif(n, 0, 30)
d$status <- as.numeric(d$time <= cc); d$time <- pmin(d$time, cc); d$w <- runif(n, 0.3, 3)

# S1: MB weights at / beyond the last time point
tp <- 1:10; S <- seq(0.95, 0.5, length.out = 10)
w_eq <- wt.rg.S(S, "MB", mb_tstar = 10, tpoints = tp); w_beyond <- wt.rg.S(c(S[-10], 0), "MB", mb_tstar = 11, tpoints = tp)
say("S1", "mb_tstar == max(tpoints): unique weights = ", paste(unique(w_eq), collapse = ","),
    " | mb_tstar > max(tpoints) with S reaching 0: weights = ", paste(format(utils::tail(w_beyond, 3), digits = 4), collapse = ", "),
    " (last S_left = ", utils::tail(c(1, S[-10]), 1), ")")
S0 <- c(S[-10], 0); w_b2 <- wt.rg.S(c(S0[-10], 0, 0)[1:10], "MB", mb_tstar = 11, tpoints = tp)
say("S1", "S_left = 0 inside the grid, mb_tstar beyond: any infinite weight = ", any(is.infinite(wt.rg.S(c(0.9, 0.5, 0, 0, 0), "MB", mb_tstar = 99, tpoints = 1:5))))

# S2: requires_cols typo -- what happens when survG is absent with scheme = "schemper"
dc <- df_counting(d, "time", "status", "treat", check.KM = FALSE, get.cox = FALSE)
dc2 <- dc; dc2$survG <- NULL
r <- run(wlr_dhat_estimates(dc2, scheme = "schemper", tzero = 10))
say("S2", "wlr_dhat_estimates(scheme='schemper') without survG -> ", if (is.null(r$error)) "no error" else paste("error:", r$error))
r <- run(wlr_cumulative(dc2[c("at_points", "nbar0", "nbar1", "ybar0", "ybar1", "survP")], scheme = "schemper"))
say("S2", "note: dfcounting$survG partially matches dfcounting$survG_all (length ", length(dc2$survG), ") once survG is removed")
say("S2", "wlr_cumulative(scheme='schemper') without survG -> ", if (is.null(r$error)) "no error" else paste("error:", r$error))

# S3: default at_points uses raw times, estimation uses aeqSurv-fixed times
dn <- d; dn$time[2] <- dn$time[1] + 1e-12 * dn$time[1]; dn$status[1:2] <- 1; dn$treat[1:2] <- 0
fixed <- survival::aeqSurv(survival::Surv(dn$time, dn$status))[, "time"]
k <- KM_diff(dn, "time", "status", "treat")
say("S3", "near-tie pair raw times differ by ", signif(diff(dn$time[1:2]), 3), "; fixed times differ by ", diff(fixed[1:2]),
    "; KM_diff at_points contains both raw values: ", all(dn$time[1:2] %in% k$at_points))
ap <- k$at_points; i <- which(ap %in% dn$time[1:2])
say("S3", "S0 at the two raw points: ", paste(format(k$surv0[i], digits = 12), collapse = " vs "))

# S7: get_validated_weights writes with cat() when mb_tstar is missing
r <- run(get_validated_weights(data.frame(at_points = 1:5, S.pool = seq(1, .6, length.out = 5)), scheme = "MB", scheme_params = list()))
say("S7", "stdout: '", paste(r$stdout, collapse = " "), "'; then error: ", r$error)

# S8: strata.name not validated; event coding not validated in df_counting
r <- run(df_counting(d, "time", "status", "treat", strata.name = "no_such_column"))
say("S8", "df_counting(strata.name = 'no_such_column') -> ", if (is.null(r$error)) "no error" else paste("error:", r$error))
d3 <- d; d3$status[1] <- 3
r <- run(df_counting(d3, "time", "status", "treat", get.cox = FALSE))
say("S8", "df_counting with an event coded 3 -> ", if (is.null(r$error)) paste("no error;", length(r$warnings), "warnings") else paste("error:", r$error))

# S9: same seedstart for both arms / every stratum
g0 <- d[d$treat == 0, ]; g1 <- d[d$treat == 1, ]
set.seed(5); a <- calculate_risk_event_counts(g0$time, g0$status, g0$w, sort(g0$time), draws = 50, seedstart = 1)
set.seed(6); b <- calculate_risk_event_counts(g0$time, g0$status, g0$w, sort(g0$time), draws = 50, seedstart = 1)
say("S9", "calculate_risk_event_counts(draws=50) identical across calls regardless of caller seed: ", identical(a, b),
    " (set.seed(seedstart) at R/df_counting_helpers.R:261; df_counting passes the same seedstart to every arm/stratum call)")

# N2: df_counting unweighted, z.score NaN -> 'missing value where TRUE/FALSE needed'
d1 <- data.frame(time = c(3, 5), status = c(0, 1), treat = c(0, 1))
msg <- where(df_counting(d1, "time", "status", "treat", check.KM = FALSE, get.cox = FALSE))
say("N2", "n = 1 per arm, unweighted: ", msg, " | call stack tail: ", paste(.last_calls, collapse = " <- "))
dc1 <- run(df_counting(transform(d, status = 0), "time", "status", "treat", check.KM = FALSE, get.cox = FALSE))
say("N2", "all censored, check.KM = FALSE: ", if (is.null(dc1$error)) "no error" else dc1$error)

# N3: singular system from survdiff on degenerate data
dh <- data.frame(time = c(8, 8), status = c(1, 1), treat = c(0, 1))
msg <- where(df_counting(dh, "time", "status", "treat", check.KM = FALSE, get.cox = FALSE))
say("N3", "n = 1 per arm, tied event times: ", msg, " | call stack tail: ", paste(.last_calls, collapse = " <- "))

# N5: cox_rhogamma when root finding fails
dz <- d; dz$status[dz$treat == 1] <- 0
dcz <- df_counting(dz, "time", "status", "treat", check.KM = FALSE, get.cox = FALSE)
r <- run(cox_rhogamma(dcz, scheme = "fh", scheme_params = list(rho = 0, gamma = 0)))
say("N5", "cox_rhogamma with no events in arm 1 -> ", if (is.null(r$error)) paste("no error; bhat =", r$value$fit$bhat) else paste("error:", r$error))

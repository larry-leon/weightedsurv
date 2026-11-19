#!/usr/bin/env Rscript
# n_findings_after.R -- after-state of Stage-1 new findings N1-N9 (Stage 2, Step 5).
# Run from the repo root: Rscript dev/audit/stage2/n_findings_after.R > dev/audit/stage2/n_findings_after.txt 2>&1
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
grDevices::pdf(NULL)
run <- function(expr) {
  ws <- character(0)
  val <- tryCatch(withCallingHandlers(expr, warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }),
                  error = function(e) structure(conditionMessage(e), class = "err"))
  list(value = if (inherits(val, "err")) NULL else val, error = if (inherits(val, "err")) unclass(val) else NULL, warnings = ws)
}
say <- function(id, r, extra = "") cat(sprintf("[%s] %s | warnings: %d%s%s\n", id,
  if (is.null(r$error)) "result returned" else paste("error:", r$error), length(r$warnings),
  if (length(r$warnings)) paste0(" (", paste(unique(r$warnings), collapse = " / "), ")") else "", extra))
set.seed(1)
n <- 120; d <- data.frame(treat = rep(0:1, each = n / 2)); d$time <- rexp(n, 0.08); cc <- runif(n, 0, 30)
d$status <- as.numeric(d$time <= cc); d$time <- pmin(d$time, cc); d$w <- runif(n, 0.3, 3)

r <- run(df_counting(transform(d, status = 0), "time", "status", "treat")); say("N1", r, " -- all censored, defaults (check.KM = TRUE, get.cox = TRUE)")
r <- run(df_counting(data.frame(time = c(3, 5), status = c(0, 1), treat = c(0, 1)), "time", "status", "treat", check.KM = FALSE, get.cox = FALSE))
say("N2", r, sprintf(" -- n = 1 per arm, unweighted; z.score = %s", if (is.null(r$value)) "-" else format(r$value$z.score)))
r <- run(df_counting(data.frame(time = c(8, 8), status = c(1, 1), treat = c(0, 1)), "time", "status", "treat", check.KM = FALSE, get.cox = FALSE))
say("N3", r, sprintf(" -- n = 1 per arm, tied event times; logrank chisq = %s", if (is.null(r$value)) "-" else format(r$value$logrank_results$chisq)))
dz <- d; dz$status[dz$treat == 1] <- 0
dcz <- df_counting(dz, "time", "status", "treat", check.KM = FALSE, get.cox = FALSE)
r <- run(cox_rhogamma(dcz, scheme = "fh", scheme_params = list(rho = 0, gamma = 0)))
say("N5", r, sprintf(" -- no events in arm 1; bhat = %s", if (is.null(r$value)) "-" else format(r$value$fit$bhat)))
dw <- d; dw$w[1:10] <- 0
r <- run(df_counting(dw, "time", "status", "treat", weight.name = "w")); say("N6", r, sprintf(" -- zero weights, get.cox = TRUE; HR = %s", if (is.null(r$value)) "-" else format(unname(r$value$cox_results$hr))))
dn <- d; dn$status[dn$treat == 0] <- 0
r <- run(plotKM.band_subgroups(dn, "time", "status", "treat", draws.band = 0)); say("N7", r, " -- arm without events, Maxtau = NULL")
r <- run(plotKM.band_subgroups(dn, "time", "status", "treat", draws.band = 0, Maxtau = 20)); say("N7", r, " -- same, Maxtau = 20")
r <- run(df_counting(d, "time", "status", "treat", weight.name = "w"))
u <- survival::survdiff(survival::Surv(time, status) ~ treat, data = d)$chisq
cat(sprintf("[N8] weighted df_counting: logrank_results$chisq = %.6f; unweighted survdiff chisq = %.6f (unchanged: survdiff takes no case weights)\n", r$value$logrank_results$chisq, u))
dneg <- d; dneg$time[1] <- -1
for (f in c("KM_diff", "df_counting", "plotKM.band_subgroups")) {
  r <- run(switch(f, KM_diff = KM_diff(dneg, "time", "status", "treat"), df_counting = df_counting(dneg, "time", "status", "treat"),
                  plotKM.band_subgroups = plotKM.band_subgroups(dneg, "time", "status", "treat", draws.band = 0)))
  say("N9", r, paste(" --", f, "with a negative time"))
}

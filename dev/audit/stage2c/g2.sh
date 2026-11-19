S=dev/audit/stage2c

# 1. branch; master untouched
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$MASTER_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"

# 2. API identical; vignettes untouched
Rscript $S/api.R compare && echo "OK 2a" || echo "FAIL: API changed"
git diff --quiet "$START_SHA" HEAD -- vignettes/ && echo "OK 2b" || echo "FAIL: vignettes changed"

# 3. golden: identical outside the Stage-2c whitelist; contract checks pass
Rscript $S/golden.R compare final && echo "OK 3" || echo "FAIL: golden"

# 4. R1: weighted SE equals survfit robust (ties included), is scale-free, and is Greenwood at equal weights
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
  set.seed(1); n <- 150; U <- ceiling(rexp(n, .1) * 2) / 2; D <- rbinom(n, 1, .7); w <- runif(n, .3, 3); at <- sort(unique(U))
  se <- function(w) { re <- calculate_risk_event_counts(U, D, w, at); sqrt(pmax(KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)$sig2_KM, 0)) }
  y  <- survival::Surv(U, D)
  rb <- summary(survival::survfit(y ~ 1, weights = w, robust = TRUE), times = at, extend = TRUE)$std.err
  gw <- summary(survival::survfit(y ~ 1), times = at, extend = TRUE)$std.err
  k  <- is.finite(rb) & rb > 0 & is.finite(gw) & gw > 0
  stopifnot(max(abs(se(w)[k] / rb[k] - 1)) < 1e-8, max(abs(se(w / sum(w))[k] / rb[k] - 1)) < 1e-8, max(abs(se(10 * w)[k] / rb[k] - 1)) < 1e-8,
            max(abs(se(rep(1, n))[k] / gw[k] - 1)) < 1e-8, max(abs(se(rep(2, n))[k] / gw[k] - 1)) < 1e-8); cat("OK 4\n")'

# 5. R2: variance of the resampled weighted curves agrees with the analytic variance
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
  set.seed(2); N <- 300; x <- rnorm(N); e <- plogis(-0.5 + x); a <- rbinom(N, 1, e)
  T <- rexp(N, 0.08 * exp(0.7 * x - 0.4 * a)); cc <- runif(N, 0, 40)
  d <- data.frame(time = pmin(T, cc), status = as.numeric(T <= cc), treat = a, w = ifelse(a == 1, 1 / e, 1 / (1 - e)))
  f <- KM_diff(d, "time", "status", "treat", weight.name = "w", draws.band = 4000, show_resamples = FALSE)
  for (g in 0:1) { S <- f[[paste0("surv", g)]]; v <- f[[paste0("sig2_surv", g)]]; vs <- apply(f[[paste0("surv", g, "_star")]], 1, var)
    k <- S > 0.1 & S < 0.9 & v > 0; r <- vs[k] / v[k]; stopifnot(sum(k) > 20, min(r) > 0.85, max(r) < 1.15) }
  cat("OK 5\n")'

# 6. seed probe: only the default RMST rule still differs from survfit
Rscript dev/audit/km_audit_probe.R --no-write > $S/probe_2c.txt 2>&1
[ "$(awk -F'|' '/^\|/ && $3 ~ /^ REPRODUCED/ {gsub(/ /,"",$2); print $2}' $S/probe_2c.txt | sort | tr '\n' ' ')" = "I3 " ] && echo "OK 6a" || echo "FAIL: probe statuses"
[ "$(grep -c '| PASS ' $S/probe_2c.txt)" = "3" ] && echo "OK 6b" || echo "FAIL: controls"

# 7. tests and check
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)' > $S/tests_2c.txt 2>&1 && echo "OK 7a" || echo "FAIL: tests"
grep -E '^Status: ' $S/00check_2c.log | grep -Evq 'ERROR|WARNING' && echo "OK 7b" || echo "FAIL: R CMD check"

# 8. log-rank and Cox code untouched
git diff --quiet "$START_SHA" HEAD -- R/weightedCox.R && echo "OK 8" || echo "FAIL: weightedCox.R changed"

# 9. deliverables and record
for f in golden_s2b.rds api_s2b.rds probe_2c.txt 00check_2c.log; do [ -f "$S/$f" ] || echo "FAIL: missing $S/$f"; done
[ -f dev/audit/REPORT_km_fixes_stage2c.md ] && echo "OK 9 (if no FAIL above)"
git log --oneline "$START_SHA"..HEAD

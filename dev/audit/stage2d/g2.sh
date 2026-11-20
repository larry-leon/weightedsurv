S=dev/audit/stage2d

# 1. branch; master untouched
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$MASTER_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"

# 2. API: the B12 default and nothing else; vignettes untouched
Rscript $S/api.R compare && echo "OK 2a" || echo "FAIL: API"
git diff --quiet "$START_SHA" HEAD -- vignettes/ && echo "OK 2b" || echo "FAIL: vignettes changed"

# 3. golden: identical outside the Stage-2d whitelist; coxph checks pass
Rscript $S/golden.R compare final && echo "OK 3" || echo "FAIL: golden"

# 4. L1, L2: against coxph (ties and strata included), scale-free, unchanged at equal weights
Rscript -e 'pkgload::load_all(".", quiet = TRUE); suppressMessages(library(survival))
  set.seed(4); n <- 300; x <- rnorm(n); e <- plogis(-0.5 + x); a <- rbinom(n, 1, e)
  T <- rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- runif(n, 0, 60)
  mk <- function(ties) data.frame(time = if (ties) ceiling(pmin(T, cc)) else pmin(T, cc), status = as.numeric(T <= cc), treat = a,
                                  w = ifelse(a == 1, 1 / e, 1 / (1 - e)), s = ifelse(x > 0, "hi", "lo"), one = 1)
  q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
  for (ties in c(FALSE, TRUE)) { d <- mk(ties); d$w10 <- 10 * d$w; d$wn <- d$w / sum(d$w)
    dc  <- q(df_counting(d, "time", "status", "treat", weight.name = "w"))
    cx  <- q(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 0))
    dcs <- q(df_counting(d, "time", "status", "treat", weight.name = "w", strata.name = "s"))
    cp  <- coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow")
    cps <- coxph(Surv(time, status) ~ treat + strata(s), data = d, weights = w, robust = TRUE, ties = "breslow")
    z10 <- q(df_counting(d, "time", "status", "treat", weight.name = "w10"))$z.score
    zn  <- q(df_counting(d, "time", "status", "treat", weight.name = "wn"))$z.score
    u   <- q(df_counting(d, "time", "status", "treat")); c1 <- q(df_counting(d, "time", "status", "treat", weight.name = "one"))
    stopifnot(abs(dc$z.score^2 / cp$rscore - 1) < 1e-8, abs(cx$fit$z.score^2 / cp$rscore - 1) < 1e-8,
              abs(cx$fit$sig_bhat_asy / sqrt(cp$var) - 1) < 1e-6, abs(dcs$z.score_stratified^2 / cps$rscore - 1) < 1e-8,
              is.finite(zn), abs(z10 / dc$z.score - 1) < 1e-9, abs(zn / dc$z.score - 1) < 1e-9,
              identical(u$z.score, c1$z.score), identical(u$sig2_lr, c1$sig2_lr)) }
  cat("OK 4\n")'

# 5. L3: the resampling SE agrees with the analytic SE
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
  set.seed(5); n <- 300; x <- rnorm(n); e <- plogis(-0.5 + x); a <- rbinom(n, 1, e)
  T <- rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- runif(n, 0, 60)
  d <- data.frame(time = pmin(T, cc), status = as.numeric(T <= cc), treat = a, w = ifelse(a == 1, 1 / e, 1 / (1 - e)))
  q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
  dc <- q(df_counting(d, "time", "status", "treat", weight.name = "w"))
  cx <- q(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 4000))
  r <- cx$fit$sig_bhat_star / cx$fit$sig_bhat_asy; stopifnot(is.finite(r), r > 0.95, r < 1.05); cat("OK 5\n")'

# 6. seed probe: statuses as after Stage 2c
Rscript dev/audit/km_audit_probe.R --no-write > $S/probe_2d.txt 2>&1
[ "$(awk -F'|' '/^\|/ && $3 ~ /^ REPRODUCED/ {gsub(/ /,"",$2); print $2}' $S/probe_2d.txt | sort | tr '\n' ' ')" = "I3 " ] && echo "OK 6a" || echo "FAIL: probe statuses"
[ "$(grep -c '| PASS ' $S/probe_2d.txt)" = "3" ] && echo "OK 6b" || echo "FAIL: controls"

# 7. tests and check
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)' > $S/tests_2d.txt 2>&1 && echo "OK 7a" || echo "FAIL: tests"
grep -E '^Status: ' $S/00check_2d.log | grep -Evq 'ERROR|WARNING' && echo "OK 7b" || echo "FAIL: R CMD check"

# 8. deliverables and record
for f in golden_s2c.rds api_s2c.rds probe_2d.txt 00check_2d.log; do [ -f "$S/$f" ] || echo "FAIL: missing $S/$f"; done
[ -f dev/audit/REPORT_km_fixes_stage2d.md ] && echo "OK 8 (if no FAIL above)"
git log --oneline "$START_SHA"..HEAD

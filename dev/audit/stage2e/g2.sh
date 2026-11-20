#!/usr/bin/env bash
# g2.sh -- Gate G2 of TASK_km_fixes_stage2e.md (section 8), the block as given in the task document.
# Run from the repo root with START_SHA and MASTER_SHA set:   bash dev/audit/stage2e/g2.sh   (checks 1-12)
#                                                             bash dev/audit/stage2e/g2.sh 13 (check 13, after the final commit)
if [ "$1" = "13" ]; then
# 13. after the final commit
[ -z "$(git status --porcelain)" ] && echo "OK 13" || echo "FAIL: uncommitted changes"
exit 0
fi
S=dev/audit/stage2e

# 1. branch; master untouched
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$MASTER_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"

# 2. API: the listed differences and nothing else; vignettes and README untouched
Rscript $S/api.R compare && echo "OK 2a" || echo "FAIL: API"
git diff --quiet "$START_SHA" HEAD -- vignettes/ README.md && echo "OK 2b" || echo "FAIL: vignettes or README changed"

# 3. golden: identical outside the Stage-2e whitelist; new elements only those listed
Rscript $S/golden.R compare final && echo "OK 3" || echo "FAIL: golden"

# 4. Rd: zero-weight rows are inert
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(survival))
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
set.seed(3); n <- 300; x <- rnorm(n); a <- rbinom(n, 1, plogis(-0.3 + x))
T <- rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- runif(n, 0, 60); keep <- runif(n) < 0.7
for (ties in c(FALSE, TRUE)) {
  d <- data.frame(time = if (ties) ceiling(pmin(T, cc)) else pmin(T, cc), status = as.numeric(T <= cc), treat = a, s = ifelse(x > 0, "hi", "lo"))
  d$w01 <- as.numeric(keep); d$wu <- ifelse(keep, ifelse(a == 1, 1 / plogis(-0.3 + x), 1 / (1 - plogis(-0.3 + x))), 0)
  sub <- d[keep, ]
  for (wn in c("w01", "wu")) {
    f <- q(df_counting(d, "time", "status", "treat", weight.name = wn, strata.name = "s"))
    p <- q(df_counting(sub, "time", "status", "treat", weight.name = if (wn == "w01") NULL else wn, strata.name = "s"))
    cf <- q(cox_rhogamma(f, scheme = "fh", scheme_params = list(rho = 0, gamma = 0.5)))
    cp <- q(cox_rhogamma(p, scheme = "fh", scheme_params = list(rho = 0, gamma = 0.5)))
    kf <- q(KM_diff(d, "time", "status", "treat", weight.name = wn, at_points = p$at_points))
    kp <- q(KM_diff(sub, "time", "status", "treat", weight.name = if (wn == "w01") NULL else wn, at_points = p$at_points))
    j <- match(p$at_points, f$at_points)
    stopifnot(!anyNA(j), max(abs(f$sig2_surv1[j] - p$sig2_surv1)) < 1e-12, max(abs(f$sig2_surv0[j] - p$sig2_surv0)) < 1e-12,
              abs(f$z.score - p$z.score) < 1e-10, abs(f$z.score_stratified - p$z.score_stratified) < 1e-10,
              abs(cf$fit$sig_bhat_asy / cp$fit$sig_bhat_asy - 1) < 1e-10, abs(cf$z.score - cp$z.score) < 1e-10,
              max(abs(kf$sig2_dhat - kp$sig2_dhat)) < 1e-12)
  }
}
cat("OK 4\n")'

# 5. Rb, Ra: the survival-based outputs with weights; effective sample size
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(survival))
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
set.seed(4); n <- 300; x <- rnorm(n); e <- plogis(-0.5 + x); a <- rbinom(n, 1, e)
T <- rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- runif(n, 0, 60)
ess <- function(w) sum(w)^2 / sum(w^2)
for (ties in c(FALSE, TRUE)) {
  d <- data.frame(time = if (ties) ceiling(pmin(T, cc)) else pmin(T, cc), status = as.numeric(T <= cc), treat = a,
                  w = ifelse(a == 1, 1 / e, 1 / (1 - e)), s = ifelse(x > 0, "hi", "lo"), two = 2)
  dc  <- q(df_counting(d, "time", "status", "treat", weight.name = "w"))
  dcs <- q(df_counting(d, "time", "status", "treat", weight.name = "w", strata.name = "s"))
  u   <- q(df_counting(d, "time", "status", "treat")); c2 <- q(df_counting(d, "time", "status", "treat", weight.name = "two"))
  cpb <- coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow")
  cps <- coxph(Surv(time, status) ~ treat + strata(s), data = d, weights = w, robust = TRUE, ties = "breslow")
  cpe <- summary(coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE))
  stopifnot(abs(dc$logrank_results$chisq / cpb$rscore - 1) < 1e-8, abs(dc$logrank_results$chisq / dc$z.score^2 - 1) < 1e-8,
            abs(dcs$logrank_results$chisq / cps$rscore - 1) < 1e-8, abs(dcs$logrank_results$chisq / dcs$z.score_stratified^2 - 1) < 1e-8,
            abs(dc$cox_results$score / cpe$robscore[1] - 1) < 1e-8,
            abs(u$logrank_results$chisq / survdiff(Surv(time, status) ~ treat, data = d)$chisq - 1) < 1e-10,
            identical(u$logrank_results$chisq, c2$logrank_results$chisq),
            abs(dc$ess$ess[1] / ess(d$w[d$treat == 1]) - 1) < 1e-10, abs(dc$ess$ess[2] / ess(d$w[d$treat == 0]) - 1) < 1e-10,
            all(abs(u$ess$ess - u$ess$n) < 1e-9), all(abs(c2$ess$ess - c2$ess$n) < 1e-9), all(abs(dc$ess$ratio - dc$ess$ess / dc$ess$n) < 1e-12))
}
cat("OK 5\n")'

# 6. W2: Kaplan-Meier with estimated weights against the n x m reference; weight form; reduction to known; misuse; resampling
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(survival))
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
err <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
set.seed(11); n <- 400; x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4); A <- rbinom(n, 1, plogis(-0.4 + 0.9 * x1 - 0.5 * x2))
T <- rexp(n, 0.04 * exp(0.5 * x1 - 0.4 * A)); C <- runif(n, 0, 72)
# n_g x m matrix of Kaplan-Meier influence values of one group on a sorted grid holding its event times
kmpsi <- function(U, D, W, g) {
  ybar <- colSums(outer(U, g, ">=") * W); dN <- diff(c(0, colSums(outer(U[D == 1], g, "<=") * W[D == 1])))
  ok <- ybar > 0 & (ybar - dN) > 1e-10 * ybar
  t(apply(sweep(outer(U, g, "==") * D, 2, ifelse(ok, 1 / (ybar - dN), 0), "*") - sweep(outer(U, g, ">="), 2, ifelse(ok, dN / (ybar * (ybar - dN)), 0), "*"), 1, cumsum))
}
for (ties in c(FALSE, TRUE)) {
  d <- data.frame(time = if (ties) ceiling(pmin(T, C)) else pmin(T, C), status = as.numeric(T <= C), treat = A, x1 = x1, x2 = x2)
  pw <- ps_weights(treat ~ x1 + x2, data = d, type = "ipw"); pst <- ps_weights(treat ~ x1 + x2, data = d, type = "stabilized")
  d$w <- pw$weights; d$wst <- pst$weights; g <- sort(unique(d$time))
  fk <- q(KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = g, draws.band = 300))
  fe <- q(KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = g, draws.band = 300, weight_variance = "estimated", weight_model = pw))
  fs <- q(KM_diff(d, "time", "status", "treat", weight.name = "wst", at_points = g, draws.band = 300, weight_variance = "estimated", weight_model = pst))
  pz <- pw; pz$IF <- pw$IF * 0
  f0 <- q(KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = g, draws.band = 300, weight_variance = "estimated", weight_model = pz))
  de <- q(df_counting(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw)); jd <- match(g, de$at_points)
  # reference: Phi_a[i, j] = 1(i in arm a) w_i psi_i(g_j) + IF_i (sum_{k in a} psi_k(g_j) dw_k)
  Phi <- list()
  for (a in c(0, 1)) { ia <- which(d$treat == a); P <- kmpsi(d$time[ia], d$status[ia], d$w[ia], g)
    Fa <- pw$IF %*% t(crossprod(P, pw$dw[ia, , drop = FALSE])); Fa[ia, ] <- Fa[ia, ] + d$w[ia] * P; Phi[[a + 1]] <- Fa }
  v0 <- fe$surv0^2 * colSums(Phi[[1]]^2); v1 <- fe$surv1^2 * colSums(Phi[[2]]^2)
  vd <- colSums((Phi[[2]] * rep(fe$surv1, each = n) - Phi[[1]] * rep(fe$surv0, each = n))^2)
  stopifnot(identical(fe$dhat, fk$dhat), identical(fe$surv0, fk$surv0), identical(fe$surv1, fk$surv1),
            rel(fe$sig2_surv0, v0) < 1e-8, rel(fe$sig2_surv1, v1) < 1e-8, rel(fe$sig2_dhat, vd) < 1e-8,
            rel(de$sig2_surv0[jd], v0) < 1e-8, rel(de$sig2_surv1[jd], v1) < 1e-8,
            rel(fs$sig2_dhat, fe$sig2_dhat) < 1e-8, rel(fs$dhat_star, fe$dhat_star) < 1e-8, abs(fs$c_alpha_band / fe$c_alpha_band - 1) < 1e-8,
            rel(f0$sig2_dhat, fk$sig2_dhat) < 1e-10, rel(f0$dhat_star, fk$dhat_star) < 1e-10,
            identical(fe$weight_variance, "estimated"), identical(fk$weight_variance, "known"), identical(de$weight_variance, "estimated"),
            err(KM_diff(d, "time", "status", "treat", weight_variance = "estimated", weight_model = pw)),
            err(KM_diff(d, "time", "status", "treat", weight.name = "wst", weight_variance = "estimated", weight_model = pw)),
            err(KM_diff(d[-1, ], "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw)),
            err(KM_diff(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw$model)),
            err(KM_diff(d, "time", "status", "treat", weight.name = "w", draws = 50, weight_variance = "estimated", weight_model = pw)),
            err(df_counting(d, "time", "status", "treat", weight.name = "wst", weight_variance = "estimated", weight_model = pw)))
  if (!ties) {
    f4 <- q(KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = g, draws.band = 4000, weight_variance = "estimated", weight_model = pw))
    pos <- f4$sig2_dhat > 0; r <- median(apply(f4$dhat_star[pos, ], 1, var) / f4$sig2_dhat[pos])
    rm <- q(cumulative_rmst_bands(d, f4, "time", "status", "treat", weight.name = "w", draws_sb = 300, plot = FALSE, weight_variance = "estimated", weight_model = pw))
    stopifnot(r > 0.95, r < 1.05, identical(rm$weight_variance, "estimated"),
              err(cumulative_rmst_bands(d, fk, "time", "status", "treat", weight.name = "w", draws_sb = 300, plot = FALSE, weight_variance = "estimated", weight_model = pw)))
    cat(sprintf("resampled / analytic variance of the difference (median over the grid): %.3f\n", r))
  }
}
cat("OK 6\n")'

# 7. W3, W4: log-rank and Cox against the coxph residual routes (ties, strata); cumulative; resampling
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(survival))
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
set.seed(11); n <- 400; x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4); A <- rbinom(n, 1, plogis(-0.4 + 0.9 * x1 - 0.5 * x2))
T <- rexp(n, 0.04 * exp(0.5 * x1 - 0.4 * A)); C <- runif(n, 0, 72)
P00 <- list(rho = 0, gamma = 0)
for (ties in c(FALSE, TRUE)) for (type in c("ipw", "stabilized")) {
  d <- data.frame(time = if (ties) ceiling(pmin(T, C)) else pmin(T, C), status = as.numeric(T <= C), treat = A, x1 = x1, x2 = x2, s = ifelse(x2 == 1, "hi", "lo"))
  pw <- ps_weights(treat ~ x1 + x2, data = d, type = type); d$w <- pw$weights
  adj <- function(r) { r <- as.numeric(r); as.numeric(r + pw$IF %*% colSums((r / d$w) * pw$dw)) }
  sres <- function(f) { r <- as.numeric(residuals(f, type = "score", weighted = TRUE)); sum(r)^2 / sum(adj(r)^2) }
  chi  <- sres(coxph(Surv(time, status) ~ treat, data = d, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow"))
  chis <- sres(coxph(Surv(time, status) ~ treat + strata(s), data = d, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow"))
  se_b <- sqrt(sum(adj(residuals(coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow"), type = "dfbeta", weighted = TRUE))^2))
  se_e <- sqrt(sum(adj(residuals(coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE), type = "dfbeta", weighted = TRUE))^2))
  dk  <- q(df_counting(d, "time", "status", "treat", weight.name = "w"))
  de  <- q(df_counting(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw))
  des <- q(df_counting(d, "time", "status", "treat", weight.name = "w", strata.name = "s", weight_variance = "estimated", weight_model = pw))
  ck  <- q(cox_rhogamma(dk, scheme = "fh", scheme_params = P00)); cx <- q(cox_rhogamma(de, scheme = "fh", scheme_params = P00))
  one <- q(wlr_cumulative(de, scheme = "fh", scheme_params = P00)); cum <- q(wlr_cumulative(de, scheme = "fh", scheme_params = P00, return_cumulative = TRUE))
  t0 <- median(d$time); j0 <- max(which(de$at_points <= t0)); dcen <- d; dcen$status <- d$status * (d$time <= t0); dcen$time <- pmin(d$time, t0)
  dn <- q(df_counting(dcen, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw))
  stopifnot(identical(de$lr, dk$lr), identical(cx$fit$bhat, ck$fit$bhat), identical(de$cox_results$hr, dk$cox_results$hr),
            abs(de$z.score^2 / chi - 1) < 1e-8, abs(des$z.score_stratified^2 / chis - 1) < 1e-8,
            abs(de$logrank_results$chisq / chi - 1) < 1e-8, abs(des$logrank_results$chisq / chis - 1) < 1e-8,
            abs(cx$fit$z.score^2 / chi - 1) < 1e-8, abs(cx$fit$sig_bhat_asy / se_b - 1) < 1e-6, abs(de$cox_results$se / se_e - 1) < 1e-8,
            abs(one$sig2.score / de$sig2_lr - 1) < 1e-10, abs(cum$sig2.score[length(cum$sig2.score)] / de$sig2_lr - 1) < 1e-10,
            abs(cum$sig2.score[j0] / dn$sig2_lr - 1) < 1e-8,
            identical(de$weight_variance, "estimated"), identical(cx$weight_variance, "estimated"), identical(ck$weight_variance, "known"))
  if (!ties && type == "ipw") {
    cs <- q(cox_rhogamma(de, scheme = "fh", scheme_params = P00, draws = 4000)); r <- cs$fit$sig_bhat_star / cs$fit$sig_bhat_asy
    stopifnot(is.finite(r), r > 0.95, r < 1.05, cx$fit$sig_bhat_asy < ck$fit$sig_bhat_asy)
    cat(sprintf("Cox: resampling SE / analytic SE %.3f; analytic SE estimated %.4f, known %.4f\n", r, cx$fit$sig_bhat_asy, ck$fit$sig_bhat_asy))
  }
}
cat("OK 7\n")'

# 8. W1: ps_weights() against numerical derivatives
Rscript -e 'pkgload::load_all(".", quiet = TRUE)
err <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
set.seed(11); n <- 400; x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4); A <- rbinom(n, 1, plogis(-0.4 + 0.9 * x1 - 0.5 * x2))
d <- data.frame(treat = A, x1 = x1, x2 = x2)
g <- glm(treat ~ x1 + x2, family = binomial(), data = d); e <- as.numeric(fitted(g)); p1 <- mean(A); xm <- model.matrix(g)
pw <- ps_weights(treat ~ x1 + x2, data = d, type = "ipw"); pst <- ps_weights(treat ~ x1 + x2, data = d, type = "stabilized"); pm <- ps_weights(model = g, type = "ipw")
# numerical derivatives, from fits iterated to convergence
ctl <- glm.control(epsilon = 1e-13, maxit = 100); b <- coef(glm(treat ~ x1 + x2, family = binomial(), data = d, control = ctl))
wof <- function(bb) { pp <- as.numeric(plogis(xm %*% bb)); ifelse(A == 1, 1 / pp, 1 / (1 - pp)) }
num <- sapply(seq_along(b), function(j) { h <- 1e-6; bp <- b; bp[j] <- b[j] + h; bm <- b; bm[j] <- b[j] - h; (wof(bp) - wof(bm)) / (2 * h) })
ii <- c(1, 7, 50, 200, 399)
ifn <- t(sapply(ii, function(i) { wt <- rep(1, n); wt[i] <- 1 + 1e-4
  (coef(glm(treat ~ x1 + x2, family = quasibinomial(), data = d, weights = wt, control = ctl)) - b) / 1e-4 }))
stopifnot(inherits(pw, "ps_weights"), all(c("weights", "type", "treat", "IF", "dw", "model") %in% names(pw)),
          max(abs(pw$weights - ifelse(A == 1, 1 / e, 1 / (1 - e)))) < 1e-10, max(abs(pst$weights - ifelse(A == 1, p1 / e, (1 - p1) / (1 - e)))) < 1e-10,
          max(abs(pm$weights - pw$weights)) < 1e-10, max(abs(pm$IF - pw$IF)) < 1e-10, max(abs(pm$dw - pw$dw)) < 1e-10,
          all(dim(pw$IF) == c(n, 3)), all(dim(pw$dw) == c(n, 3)), all(pw$treat == A),
          max(abs(pw$dw - num)) / max(abs(num)) < 1e-5, max(abs(pw$IF[ii, ] - ifn)) / max(abs(ifn)) < 1e-3,
          max(abs(pst$dw - ifelse(A == 1, p1, 1 - p1) * pw$dw)) < 1e-12, max(abs(pst$IF - pw$IF)) < 1e-12,
          err(ps_weights(model = glm(treat ~ x1, family = binomial(link = "probit"), data = d))), err(ps_weights(model = lm(treat ~ x1, data = d))),
          err(ps_weights(treat ~ x1, data = transform(d, x1 = ifelse(seq_len(n) == 3, NA, x1)))))
invisible(utils::capture.output(print(pw)))
cat("OK 8\n")'

# 9. M7 requirements
[ -f $S/m7_required.txt ] && [ "$(grep -c '^OK' $S/m7_required.txt)" = "4" ] && ! grep -q FAIL $S/m7_required.txt && echo "OK 9" || echo "FAIL: M7 requirements"

# 10. seed probe: statuses as after Stage 2c
Rscript dev/audit/km_audit_probe.R --no-write > $S/probe_2e.txt 2>&1
[ "$(awk -F'|' '/^\|/ && $3 ~ /^ REPRODUCED/ {gsub(/ /,"",$2); print $2}' $S/probe_2e.txt | sort | tr '\n' ' ')" = "I3 " ] && echo "OK 10a" || echo "FAIL: probe statuses"
[ "$(grep -c '| PASS ' $S/probe_2e.txt)" = "3" ] && echo "OK 10b" || echo "FAIL: controls"

# 11. tests and check
Rscript -e 'testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)' > $S/tests_2e.txt 2>&1 && echo "OK 11a" || echo "FAIL: tests"
grep -E '^Status: ' $S/00check_2e.log | grep -Evq 'ERROR|WARNING' && echo "OK 11b" || echo "FAIL: R CMD check"

# 12. deliverables and record
for f in golden_s2d.rds api_s2d.rds probe_2e.txt 00check_2e.log; do [ -f "$S/$f" ] || echo "FAIL: missing $S/$f"; done
[ -f dev/audit/REPORT_km_fixes_stage2e.md ] && echo "OK 12 (if no FAIL above)"
git log --oneline "$START_SHA"..HEAD

#!/usr/bin/env Rscript
# oc_boot.R -- Part E, bootstrap benchmark (TASK_km_oc_evaluation.md, section 8).
# Cells 5, 8 and 10; the first 500 replicates; 200 resamples of subjects with replacement.
# In cell 8 the propensity model is refitted in every resample; in cells 5 and 10 the true weights
# travel with the subjects. Computed with `survival` only, on the full data and on each resample:
#   delta(18): difference of weighted survfit curves at month 18
#   R(36):     difference of the weighted restricted means to 36 (rmean = 36)
#   log HR:    coef(coxph(Surv(time, status) ~ treat, weights = w))
# boot: full-data estimate +- qnorm(0.975) * SD over resamples. boot_pct: 2.5% and 97.5% quantiles.
#
#   Rscript dev/audit/oc/oc_boot.R [--cells=5,8,10] [--reps=500] [--B=200] [--cores=12] [--out=raw]

suppressPackageStartupMessages(library(survival))
.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.here, "oc_common.R"))
opt <- parse_args(list(cells = "5,8,10", reps = "500", B = "200", cores = "12", out = "raw"))
cell_ids <- eval(parse(text = paste0("c(", opt$cells, ")")))
nrep <- as.integer(opt$reps); B <- as.integer(opt$B); cores <- as.integer(opt$cores)
QN <- c("d18", "R36", "loghr")

est3 <- function(df) {
  sf <- survfit(Surv(time, status) ~ treat, data = df, weights = w)
  s <- summary(sf, times = 18, extend = TRUE, rmean = 36)
  rmn <- s$table[, grep("rmean$", colnames(s$table))[1]]
  b <- coef(coxph(Surv(time, status) ~ treat, data = df, weights = w))
  stats::setNames(c(s$surv[2] - s$surv[1], rmn[2] - rmn[1], b), QN)
}

one_rep <- function(k, r) {
  d <- make_data(k, r); df <- d$df
  full <- suppressWarnings(est3(df))
  set.seed(6e7 + 1e6 * k + r)
  n <- nrow(df)
  bs <- matrix(NA_real_, B, 3)
  for (b in seq_len(B)) {
    db <- df[sample.int(n, n, replace = TRUE), ]
    bs[b, ] <- tryCatch(suppressWarnings({
      if (cells$wt[k] == "E") {
        eh <- fitted(glm(treat ~ X, binomial, data = db))
        db$w <- ifelse(db$treat == 1, 1 / eh, 1 / (1 - eh))
      }
      est3(db)
    }), error = function(e) rep(NA_real_, 3))
  }
  ok <- colSums(is.finite(bs))
  bs[!is.finite(bs)] <- NA
  sdv <- apply(bs, 2, sd, na.rm = TRUE)
  qs <- apply(bs, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  out <- c(stats::setNames(full, paste0("est.", QN)), stats::setNames(sdv, paste0("sd.", QN)),
           stats::setNames(qs[1, ], paste0("pct_lo.", QN)), stats::setNames(qs[2, ], paste0("pct_hi.", QN)),
           stats::setNames(ok, paste0("nok.", QN)))
  attr(out, "chk") <- d$chk
  out
}

outdir <- file.path(.here, opt$out); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("oc_boot.R cells=%s reps=%d B=%d cores=%d (survival %s)  started %s\n", opt$cells, nrep, B, cores,
            as.character(utils::packageVersion("survival")), format(Sys.time())))
t_all <- proc.time()[["elapsed"]]
for (k in cell_ids) {
  t0 <- proc.time()[["elapsed"]]
  L <- parallel::mclapply(seq_len(nrep), function(r) one_rep(k, r), mc.cores = cores)
  stopifnot(all(vapply(L, is.numeric, logical(1))))
  res <- do.call(rbind, L); rownames(res) <- seq_len(nrep)
  secs <- proc.time()[["elapsed"]] - t0
  saveRDS(list(res = res, chk = vapply(L, function(x) attr(x, "chk"), character(1)), cell_id = k, B = B, secs = secs),
          file.path(outdir, sprintf("boot_cell%02d.rds", k)))
  cat(sprintf("  cell %2d %-12s %d replicates x %d resamples %7.1f s  resamples that failed: %d\n", k, cells$cell[k], nrep, B, secs,
              sum(B - res[, "nok.d18"])))
}
cat(sprintf("oc_boot.R done: wall-clock %.1f s (%.1f min)\n", proc.time()[["elapsed"]] - t_all, (proc.time()[["elapsed"]] - t_all) / 60))

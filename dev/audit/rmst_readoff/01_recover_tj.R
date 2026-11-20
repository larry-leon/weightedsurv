#!/usr/bin/env Rscript
# 01_recover_tj.R -- TASK_rmst_readoff_truth.md, section 2, step 1: the grid time t_j <= t at which each replicate's
# RMST "at month t" was read, and the truth R(t_j), for every cell with RMST rows in the article's CSV.
#
# No package call. Data are regenerated from the recorded seeds with the generators of dev/audit/stage2e/oc_common.R
# (and, for cells 41/42, a verbatim copy of gen9() of dev/audit/ps_by_time/ps9_sim.R); each regenerated data set's
# checksum is compared with the `chk` stored in the per-replicate file.
#
# Grid and read-off rule, as coded in the scripts that produced the cells:
#   dev/audit/oc/oc_sim.R:195-196, dev/audit/stage2e/oc_sim.R:182-183 (one_rep / run_option, reused by
#   ps_by_time/ps9_sim.R): fit <- KM_diff(df, ...) with the default at_points; jj <- last_le(fit$at_points, tt)
#   oc_common.R:109: last_le(grid, t) = findInterval(t, grid) (NA when 0)
#   R/km_wlr_calculations_helpers.R:567 (identical at d9b7558:496, c195953:567, f55001f:567):
#     at_points = sort(df[[tte.name]]) -- every observed time of the data frame, events and censorings, raw times
#   R/rmst_calculations.R (cumulative_rmst_bands): at_points <- fit$at_points; rmst_time on that grid, so
#     rm$rmst_time[jj] (oc_sim.R:89, stage2e/oc_sim.R:91) is the RMST at t_j = sort(df$time)[jj].
# Routes: (a) OC cells 5, 8, 10 store R.truth.<t> = R_true(fit$at_points, k)[jj] (oc/oc_sim.R:210-211);
#         (b) M7 cells 31, 32, 35, 37, 38, 39 and new cells 41, 42: t_j from the regenerated data.
# The OC cells are also regenerated (route b as a check of the rule against the stored R(t_j) of route a).
#
#   Rscript --vanilla dev/audit/rmst_readoff/01_recover_tj.R

.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
AUD <- normalizePath(file.path(.here, ".."))
source(file.path(AUD, "stage2e", "oc_common.R"))     # make_data, gen, gen_v, m7, m7_data, R_true, R_sc, last_le, tt

# gen9() and the data of cells 41, 42: copied verbatim from dev/audit/ps_by_time/ps9_sim.R (PA, PB = its defaults)
PA <- -0.44; PB <- 0.51
gen9 <- function(N, sc) {
  X <- rnorm(N)
  H <- -log(pnorm(0.6 * X + 0.8 * rnorm(N)))
  e <- plogis(PA + PB * X)
  A <- rbinom(N, 1, e)
  T <- ifelse(A == 1, LamInv(H, 1, sc), LamInv(H, 0, sc))
  C <- runif(N, 0, 72)
  data.frame(X = X, e = e, treat = A, time = pmin(T, C), status = as.numeric(T <= C))
}
data_of <- function(src, k, r) {
  if (src == "OC") return(make_data(k, r))
  if (k %in% c(41, 42)) {
    set.seed(1e6 * k + r)
    df <- gen9(400, c(`41` = "ph", `42` = "null")[[as.character(k)]])
    return(list(df = df, chk = sprintf("%.10f", sum(df$time) + sum(df$status) + sum(df$X))))
  }
  m7_data(k, r)
}
rawfile <- function(src, k) switch(src,
  OC = file.path(AUD, "oc", "raw", sprintf("proposed_cell%02d.rds", k)),
  M7 = file.path(AUD, "stage2e", "raw", sprintf("m7_cell%02d.rds", k)),
  sim = file.path(AUD, "ps_by_time", "raw", sprintf("m7_cell%02d.rds", k)))

CELLS <- data.frame(src = c("sim", "sim", rep("M7", 6), rep("OC", 3)),
                    k = c(41, 42, 31, 32, 35, 37, 38, 39, 5, 8, 10),
                    sc = c("ph", "null", "ph", "null", "ph", "null", "null", "ph", "ph", "ph", "ph"),
                    route = c(rep("b", 8), rep("a", 3)), stringsAsFactors = FALSE)
NT <- length(tt)
out <- list(); log <- character(0)
say <- function(...) { s <- sprintf(...); log <<- c(log, s); cat(s, "\n") }
say("01_recover_tj.R  %s  R %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), getRversion())
say("months: %s", paste(tt, collapse = ", "))
for (i in seq_len(nrow(CELLS))) {
  cl <- CELLS[i, ]; st <- readRDS(rawfile(cl$src, cl$k)); n <- nrow(st$res)
  tj <- matrix(NA_real_, n, NT, dimnames = list(NULL, tt)); chk_ok <- logical(n); nmin <- 0L
  for (r in seq_len(n)) {
    d <- data_of(cl$src, cl$k, r)
    chk_ok[r] <- identical(unname(d$chk), unname(st$chk[r]))
    grid <- sort(d$df$time)                            # the default at_points of KM_diff()
    jj <- last_le(grid, tt)
    tj[r, ] <- grid[jj]
  }
  # truth at t_j: the closed form of the cell's own summary (M7 / new: R_sc; OC: R_true, the evaluation's integral)
  Rtj <- if (cl$src == "OC") {
    u <- sort(unique(tj[is.finite(tj)])); Ru <- R_true(u, cl$k)
    matrix(Ru[match(tj, u)], n, NT, dimnames = list(NULL, tt))
  } else matrix(R_sc(as.vector(tj), cl$sc), n, NT, dimnames = list(NULL, tt))
  rec <- list(src = cl$src, k = cl$k, sc = cl$sc, route = cl$route, n = n, chk_ok = chk_ok, tj = tj, Rtj_regen = Rtj)
  if (cl$route == "a") {
    Rst <- st$res[, paste0("R.truth.", tt), drop = FALSE]
    rec$Rtj <- unname(Rst); dimnames(rec$Rtj) <- list(NULL, tt)
    dd <- abs(Rst - Rtj)
    say("cell %2d (%s, route a): stored R.truth columns %d x %d, NA %d; data regenerated, checksum equal in %d of %d replicates; regenerated R_true(t_j) against stored R.truth: identical %d of %d, largest |difference| %.3g",
        cl$k, cl$src, nrow(Rst), ncol(Rst), sum(is.na(Rst)), sum(chk_ok), n, sum(Rst == Rtj, na.rm = TRUE), length(Rst), max(dd, na.rm = TRUE))
  } else {
    rec$Rtj <- Rtj
    say("cell %2d (%s, route b, %s): data regenerated from the recorded seeds, checksum equal in %d of %d replicates; t_j found at every month in %d of %d replicates",
        cl$k, cl$src, cl$sc, sum(chk_ok), n, sum(rowSums(is.finite(tj)) == NT), n)
  }
  gap <- R_sc(tt, cl$sc) - colMeans(rec$Rtj)
  say("   mean of R(t) - R(t_j) by month: %s", paste(sprintf("%g", signif(gap, 4)), collapse = ", "))
  say("   mean of t - t_j by month: %s", paste(sprintf("%.4f", tt - colMeans(tj)), collapse = ", "))
  out[[sprintf("%s_%02d", cl$src, cl$k)]] <- rec
}
# cell 31 has the data of OC cell 8 (stage2e/oc_common.R m7_data: make_data(8, r))
a <- out$OC_08; b <- out$M7_31
say("cell 31 against OC cell 8 (same data): t_j identical %s; R_sc(t_j) of cell 31 against stored R.truth of cell 8: largest |difference| %.3g",
    identical(a$tj, b$tj), max(abs(a$Rtj - b$Rtj)))
G2 <- all(vapply(out, function(x) all(x$chk_ok), logical(1)))
say("G2 route recorded for every cell: a = OC 5, 8, 10; b = 41, 42, 31, 32, 35, 37, 38, 39; regenerated checksums equal in %d of %d replicates (route b cells: %d of %d): %s",
    sum(vapply(out, function(x) sum(x$chk_ok), 0)), sum(vapply(out, function(x) x$n, 0)),
    sum(vapply(out[CELLS$route == "b"], function(x) sum(x$chk_ok), 0)), sum(vapply(out[CELLS$route == "b"], function(x) x$n, 0)),
    if (G2) "PASS" else "FAIL")
saveRDS(out, file.path(.here, "tj.rds"))
writeLines(log, file.path(.here, "01_recover_tj.txt"))
if (!G2) quit(status = 1)

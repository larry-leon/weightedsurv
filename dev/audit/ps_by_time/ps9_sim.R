#!/usr/bin/env Rscript
# ps9_sim.R -- TASK_ps_weight_results_by_time.md, section 2: the gap run (cells Es9-400 and Es9-400-null).
# Sources dev/audit/stage2e/oc_common.R and reuses, without editing them, the functions of
# dev/audit/stage2e/oc_sim.R (one_rep, run_option, run_reps, the field builders and the exact package calls): the
# function definitions of that script are evaluated here, its top-level run code is not.
#
# The two new cells are added to the table `m7` in this process only:
#   k = 41  Es9-400       PH, HR = exp(-0.4)   stabilized weights from ps_weights(treat ~ X)
#   k = 42  Es9-400-null  null, HR = 1         stabilized weights from ps_weights(treat ~ X)
# Data: as gen() of oc_common.R (same draws in the same order) with e = plogis(a + b X), a and b from
# calibrate.R; replicate r of cell k seeded 1e6 * k + r, as the other M7 cells. Per replicate: both
# weight_variance options on the same data and weights, the content of M7 cells 31/32 without their
# requirement-2 / wlr_dhat_estimates() extras.
#
#   Rscript dev/audit/ps_by_time/ps9_sim.R --pkg=<path> [--reps=2000] [--cores=12] [--out=raw] [--pilot=50]

suppressPackageStartupMessages(library(survival))
.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
S2E <- normalizePath(file.path(.here, "..", "stage2e"))
source(file.path(S2E, "oc_common.R"))
opt <- parse_args(list(pkg = NA, reps = "2000", cores = "12", out = "raw", pilot = NA, a = "-0.44", b = "0.51"))
stopifnot(!is.na(opt$pkg))
suppressMessages(pkgload::load_all(opt$pkg, quiet = TRUE, export_all = TRUE, helpers = FALSE, attach_testthat = FALSE))
PKG_VERSION <- as.character(utils::packageVersion("weightedsurv"))
PKG_SHA <- tryCatch(system(sprintf("git -C %s rev-parse HEAD", shQuote(opt$pkg)), intern = TRUE), error = function(e) NA_character_)
grDevices::pdf(NULL)
Q <- qnorm(0.975); NT <- length(tt)

# ---- the function definitions of stage2e/oc_sim.R (top-level run code not evaluated) ----
ex <- parse(file.path(S2E, "oc_sim.R"), keep.source = FALSE)
is_fun_def <- function(e) is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
fdefs <- Filter(is_fun_def, as.list(ex))
for (e in fdefs) eval(e, envir = globalenv())
FUNS_FROM_OC_SIM <- vapply(fdefs, function(e) as.character(e[[2]]), character(1))

# ---- the new cells ----
PA <- as.numeric(opt$a); PB <- as.numeric(opt$b)
new <- data.frame(k = c(41, 42), cell = c("Es9-400", "Es9-400-null"), sc = c("ph", "null"), type = "stabilized",
                  reps = as.integer(opt$reps), repro = FALSE, strat = FALSE, sub = FALSE, band = TRUE, rmst = TRUE,
                  fh01 = TRUE, both = FALSE, coxph = FALSE, stringsAsFactors = FALSE)
new$loghr_truth <- ifelse(new$sc == "ph", -0.4, 0)
m7 <- rbind(m7, new[, names(m7)])

gen9 <- function(N, sc) {                     # gen() of oc_common.R with e = plogis(PA + PB * X); same draws, same order
  X <- rnorm(N)
  H <- -log(pnorm(0.6 * X + 0.8 * rnorm(N)))
  e <- plogis(PA + PB * X)
  A <- rbinom(N, 1, e)
  T <- ifelse(A == 1, LamInv(H, 1, sc), LamInv(H, 0, sc))
  C <- runif(N, 0, 72)
  data.frame(X = X, e = e, treat = A, time = pmin(T, C), status = as.numeric(T <= C))
}
m7_data_stage2e <- m7_data
m7_data <- function(k, r) {
  if (!(k %in% c(41, 42))) return(m7_data_stage2e(k, r))
  set.seed(1e6 * k + r)
  df <- gen9(400, m7$sc[m7$k == k])
  list(df = df, g = NULL, chk = sprintf("%.10f", sum(df$time) + sum(df$status) + sum(df$X)))
}

cores <- as.integer(opt$cores)
cat(sprintf("ps9_sim.R weightedsurv %s at %s; logit e = %.2f + %.2f X; functions from stage2e/oc_sim.R: %d\n",
            PKG_VERSION, substr(PKG_SHA, 1, 7), PA, PB, length(FUNS_FROM_OC_SIM)))

# ---- pilot: --pilot=<n> replicates per cell on --cores workers, projection of the full run ----
if (!is.na(opt$pilot)) {
  np <- as.integer(opt$pilot); tot <- 0
  for (k in c(41, 42)) {
    t0 <- proc.time()[["elapsed"]]
    rr <- run_reps(k, seq_len(np), cores)
    s <- proc.time()[["elapsed"]] - t0; tot <- tot + s / np * as.integer(opt$reps)
    cat(sprintf("  pilot cell %d %-13s %d replicates on %d cores: %.1f s (%.3f s per replicate wall-clock); mean ESS ratio c / t %.3f / %.3f; ps failures %d\n",
                k, m7$cell[m7$k == k], np, cores, s, s / np, mean(rr$res[, "ess.ratio0"], na.rm = TRUE),
                mean(rr$res[, "ess.ratio1"], na.rm = TRUE), sum(rr$res[, "ps.fail"] == 1)))
  }
  cat(sprintf("projected wall-clock of the full run (2 cells x %s replicates, %d cores): %.1f min\n", opt$reps, cores, tot / 60))
  quit(status = 0)
}

# ---- run (output names as M7: <out>/m7_cellNN.rds) ----
outdir <- file.path(.here, opt$out); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("run: cells 41, 42; %s replicates; %d cores; out=%s; seeds 1e6 * k + r; started %s\n", opt$reps, cores, opt$out, format(Sys.time())))
t_all <- proc.time()[["elapsed"]]
for (k in c(41, 42)) {
  t0 <- proc.time()[["elapsed"]]
  rr <- run_reps(k, seq_len(as.integer(opt$reps)), cores)
  secs <- proc.time()[["elapsed"]] - t0
  saveRDS(c(rr, list(cell_id = k, pkg_version = PKG_VERSION, pkg_sha = PKG_SHA, secs = secs, cores = cores,
                     a = PA, b = PB, seeds = sprintf("1e6 * %d + r, r = 1..%s", k, opt$reps))),
          file.path(outdir, sprintf("m7_cell%02d.rds", k)))
  failcols <- grep("\\.fail$", colnames(rr$res), value = TRUE)
  nf <- colSums(rr$res[, failcols, drop = FALSE] == 1, na.rm = TRUE)
  cat(sprintf("  cell %d %-13s %s replicates %7.1f s  failures: %s\n", k, m7$cell[m7$k == k], opt$reps, secs,
              if (any(nf > 0)) paste(sprintf("%s=%d", names(nf)[nf > 0], nf[nf > 0]), collapse = " ") else "none"))
}
cat(sprintf("ps9_sim.R done %s: wall-clock %.1f s (%.1f min)\n", format(Sys.time()), proc.time()[["elapsed"]] - t_all,
            (proc.time()[["elapsed"]] - t_all) / 60))

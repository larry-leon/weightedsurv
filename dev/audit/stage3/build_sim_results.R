## Builds vignettes/articles/sim_results/wald_vs_logrank_full.rds from the objects of the full-size
## knit of the simulations article made in the vignette revision pass
## (dev/scratch/vignette_before/sims/env.rds, written by dev/scratch/knit_vignettes.R; not committed).
## Runs no simulation. Run from the repository root:
##   Rscript dev/audit/stage3/build_sim_results.R
src <- "dev/scratch/vignette_before/sims/env.rds"
out <- "vignettes/articles/sim_results/wald_vs_logrank_full.rds"
stopifnot(file.exists("DESCRIPTION"), file.exists(src))
e <- readRDS(src)
res_out <- e$res_out
d <- res_out$results_sims
n_scen <- length(unique(res_out$get_setup$fr$Scenario))
n_total <- n_scen * res_out$number_sims

cat("source:", src, "\n  size (bytes):", file.size(src), "\n  modified:", format(file.mtime(src), "%Y-%m-%d %H:%M:%S"),
    "\n  sha256:", system2("shasum", c("-a", "256", src), stdout = TRUE), "\n")
cat("res_out$results_sims:", nrow(d), "rows x", ncol(d), "columns; rows by scenario:",
    paste(table(d$Scenario), collapse = " "), "; NA values:", sum(is.na(d)), "\n")
cat("number_sims:", res_out$number_sims, " seedstart:", res_out$seedstart, " tminutes:", res_out$tminutes, "\n")
cat("completed:", nrow(d), " failed:", n_total - nrow(d), "\n")
stopifnot(nrow(d) == 30000, n_total - nrow(d) == 0, all(table(d$Scenario) == 5000))

## the knitted output of that run, its own completion lines
md <- readLines("dev/scratch/vignette_before/sims/sims.md")
cat("lines of the knitted output of that run:\n"); cat(paste0("  ", grep("^## (Scenarios:|Completed:|Elapsed:)", md, value = TRUE)), sep = "\n")

## columns used by the tables and figures of the article
z_keep <- grep("z", colnames(d), value = TRUE)
z_keep <- z_keep[!grepl("debiased", z_keep)]
est_keep <- c("fh00_wald", "fh01_wald", "t601_wald",
              "fh00_bhat", "fh00_bhatdebiased", "fh00_sigbhat", "fh00_sigbhatstar",
              "t601_bhat", "t601_bhatdebiased", "t601_sigbhat", "t601_sigbhatstar")
keep <- c("Scenario", z_keep, est_keep)
stopifnot(all(keep %in% colnames(d)), !anyDuplicated(keep))

r2 <- res_out
r2$results_sims <- d[, keep]
stored <- list(
  res_out = r2,
  meta = list(
    scenario_names        = unique(res_out$get_setup$fr$Name),
    trials_per_scenario   = res_out$number_sims,
    seedstart             = res_out$seedstart,
    mart_draws            = 200,   # as in the run-simulations chunk of the article at the knit (file identical to ac76f82)
    num_workers           = 12,    # same source
    n_completed           = nrow(d),
    n_failed              = n_total - nrow(d),
    elapsed_minutes       = res_out$tminutes,
    package_version       = "0.2.0",
    package_commit        = "ac76f82",   # working tree at the knit; R/ is identical at ac76f82 and at the current HEAD
    run_date              = format(as.Date(file.mtime(src))),
    columns_per_trial_run = ncol(d),
    source                = "full-size knit of this article in the vignette revision pass (dev/scratch/knit_vignettes.R)"
  )
)
dir.create(dirname(out), showWarnings = FALSE)
sizes <- c(
  all_columns_xz  = { f <- tempfile(); saveRDS(res_out, f, compress = "xz"); file.size(f) },
  kept_columns_gz = { f <- tempfile(); saveRDS(stored, f, compress = "gzip"); file.size(f) })
saveRDS(stored, out, compress = "xz")
cat("\nkept columns (", length(keep), "):", paste(keep, collapse = ", "), "\n")
cat("dropped columns (", ncol(d) - length(keep), "):", paste(setdiff(colnames(d), keep), collapse = ", "), "\n")
cat("\ncompressed size, per-trial results with all", ncol(d), "columns (xz):", sizes[["all_columns_xz"]], "bytes\n")
cat("compressed size, kept columns (gzip):", sizes[["kept_columns_gz"]], "bytes\n")
cat("written:", out, "\n  size (bytes, xz):", file.size(out), " under 5 MB:", file.size(out) < 5e6, " under 5 MiB:", file.size(out) < 5 * 1024^2, "\n")
back <- readRDS(out)
cat("read back: kept columns identical to the source columns:", identical(back$res_out$results_sims, d[, keep]), "\n")
cat("meta:\n"); str(back$meta)

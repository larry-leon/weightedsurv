#!/usr/bin/env Rscript
# =============================================================================
# api.R -- API snapshot of weightedsurv exports (TASK_km_fixes_stage2d.md, Step 1)
# Copy of dev/audit/stage2c/api.R; baseline is api_s2c.rds (the branch at the start of Stage 2d).
# compare exits 0 only when the sole difference from api_s2c.rds is the B12 default:
# plotKM.band_subgroups(draws = 0) instead of 20. With "pre" as second argument (used before B12 is
# committed) it exits 0 only when there is no difference at all.
#
#   Rscript dev/audit/stage2d/api.R make          -> api_s2c.rds (refuses to overwrite)
#   Rscript dev/audit/stage2d/api.R compare       -> exit 0 only when the B12 default is the one difference
#   Rscript dev/audit/stage2d/api.R compare pre   -> exit 0 only when there is no difference
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "make"
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
S2 <- file.path("dev", "audit", "stage2d")
API <- file.path(S2, "api_s2c.rds")

snap <- function() {
  ns <- asNamespace("weightedsurv")
  ex <- sort(getNamespaceExports(ns))
  setNames(lapply(ex, function(f) {
    fn <- get(f, envir = ns)
    fm <- formals(fn)
    list(args = names(fm),
         defaults = vapply(fm, function(x) if (identical(x, quote(expr = ))) "<none>" else paste(deparse(x), collapse = " "), ""))
  }), ex)
}

if (mode == "make") {
  if (file.exists(API)) stop(API, " exists; it is the record and is not overwritten.")
  a <- snap(); saveRDS(a, API)
  cat(sprintf("api_s2c.rds: %d exports, %d arguments\n", length(a), sum(lengths(lapply(a, `[[`, "args")))))
  quit(status = 0)
}
if (mode == "compare") {
  old <- readRDS(API); new <- snap(); bad <- character(0); appended <- character(0)
  if (!identical(names(old), names(new))) {
    bad <- c(bad, paste("exports removed:", paste(setdiff(names(old), names(new)), collapse = ", ")),
             paste("exports added:", paste(setdiff(names(new), names(old)), collapse = ", ")))
  }
  for (f in intersect(names(old), names(new))) {
    o <- old[[f]]; n <- new[[f]]; k <- length(o$args)
    if (length(n$args) < k || !identical(n$args[seq_len(k)], o$args)) {
      bad <- c(bad, sprintf("%s: arguments changed: (%s) -> (%s)", f, paste(o$args, collapse = ", "), paste(n$args, collapse = ", ")))
      next
    }
    dd <- which(n$defaults[seq_len(k)] != o$defaults)
    for (i in dd) bad <- c(bad, sprintf("%s: default of '%s' changed: %s -> %s", f, o$args[i], o$defaults[i], n$defaults[i]))
    if (length(n$args) > k) appended <- c(appended, sprintf("%s: appended %s", f,
      paste(sprintf("%s = %s", n$args[-seq_len(k)], n$defaults[-seq_len(k)]), collapse = ", ")))
  }
  allowed <- if (length(args) >= 2 && args[2] == "pre") character(0) else "plotKM.band_subgroups: default of 'draws' changed: 20 -> 0"
  found <- c(appended, bad)
  cat(sprintf("API compare: %d exports before, %d now\n", length(old), length(new)))
  cat("DIFFERENCES from api_s2c.rds:\n", if (length(found)) paste0("  ", found, collapse = "\n") else "  (none)", "\n", sep = "")
  if (setequal(found, allowed) && length(found) == length(allowed)) {
    cat(if (length(allowed)) "API: the B12 default is the only difference from api_s2c.rds\n" else "API: identical to api_s2c.rds\n"); quit(status = 0)
  }
  missing_ok <- setdiff(allowed, found); extra <- setdiff(found, allowed)
  if (length(missing_ok)) cat("ALLOWED CHANGE NOT PRESENT:\n", paste0("  ", missing_ok, collapse = "\n"), "\n", sep = "")
  if (length(extra)) cat("NOT ALLOWED:\n", paste0("  ", extra, collapse = "\n"), "\n", sep = "")
  cat("API: FAIL\n"); quit(status = 1)
}
stop("usage: api.R make | compare")

#!/usr/bin/env Rscript
# =============================================================================
# api.R -- API snapshot of weightedsurv exports (TASK_km_fixes_stage2.md, Step 2b)
#
#   Rscript dev/audit/stage2/api.R make     -> api_before.rds (refuses to overwrite)
#   Rscript dev/audit/stage2/api.R compare  -> exit 0 when the export set is unchanged and, for every
#                                             export, the old arguments are an unchanged prefix (names,
#                                             order, deparsed defaults). Appended arguments are listed.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "make"
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
S2 <- file.path("dev", "audit", "stage2")
API <- file.path(S2, "api_before.rds")

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
  cat(sprintf("api_before.rds: %d exports, %d arguments\n", length(a), sum(lengths(lapply(a, `[[`, "args")))))
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
  cat(sprintf("API compare: %d exports before, %d now\n", length(old), length(new)))
  if (length(appended)) cat("APPENDED ARGUMENTS:\n", paste0("  ", appended, collapse = "\n"), "\n", sep = "")
  if (length(bad)) { cat("API DIFFERENCES:\n", paste0("  ", bad, collapse = "\n"), "\n", sep = ""); cat("API: FAIL\n"); quit(status = 1) }
  cat("API: unchanged\n"); quit(status = 0)
}
stop("usage: api.R make | compare")

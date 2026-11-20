#!/usr/bin/env Rscript
# =============================================================================
# api.R -- API snapshot of weightedsurv exports (TASK_km_fixes_stage2e.md, Step 1)
# Copy of dev/audit/stage2d/api.R; baseline is api_s2d.rds (the branch at the start of Stage 2e).
# compare exits 0 only when the differences from api_s2d.rds are exactly those of section 5 "API" for the items
# that are done (stage2e/items_done.txt, or the comma-separated list given as second argument; "none" = no item):
#   W1: export ps_weights added;
#   W2: weight_variance = c("known", "estimated"), weight_model = NULL appended to df_counting(), KM_diff(),
#       cumulative_rmst_bands(), plotKM.band_subgroups();
#   Rc: plotKM.band_subgroups() default of draws.band 20 -> 1000.
#
#   Rscript dev/audit/stage2e/api.R make             -> api_s2d.rds (refuses to overwrite)
#   Rscript dev/audit/stage2e/api.R compare [items]  -> exit 0 only when the listed differences are the only ones
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "make"
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
S2 <- file.path("dev", "audit", "stage2e")
API <- file.path(S2, "api_s2d.rds")

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
  cat(sprintf("api_s2d.rds: %d exports, %d arguments\n", length(a), sum(lengths(lapply(a, `[[`, "args")))))
  quit(status = 0)
}
if (mode == "compare") {
  old <- readRDS(API); new <- snap(); bad <- character(0); appended <- character(0)
  for (f in setdiff(names(old), names(new))) bad <- c(bad, paste("export removed:", f))
  for (f in setdiff(names(new), names(old))) bad <- c(bad, paste("export added:", f))
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
  done <- if (length(args) >= 2) trimws(strsplit(args[2], ",")[[1]]) else {
    f <- file.path(S2, "items_done.txt"); if (file.exists(f)) trimws(readLines(f)) else character(0) }
  done <- done[nzchar(done) & done != "none"]
  ap <- 'appended weight_variance = c("known", "estimated"), weight_model = NULL'
  ALLOWED <- list(W1 = "export added: ps_weights",
                  W2 = paste0(c("df_counting", "KM_diff", "cumulative_rmst_bands", "plotKM.band_subgroups"), ": ", ap),
                  Rc = "plotKM.band_subgroups: default of 'draws.band' changed: 20 -> 1000")
  allowed <- unlist(ALLOWED[intersect(done, names(ALLOWED))], use.names = FALSE)
  if (is.null(allowed)) allowed <- character(0)
  found <- c(appended, bad)
  cat(sprintf("API compare (items done: %s): %d exports before, %d now\n", if (length(done)) paste(done, collapse = ",") else "none", length(old), length(new)))
  cat("DIFFERENCES from api_s2d.rds:\n", if (length(found)) paste0("  ", found, collapse = "\n") else "  (none)", "\n", sep = "")
  if (setequal(found, allowed) && length(found) == length(allowed)) {
    cat(if (length(allowed)) "API: the listed differences are the only differences from api_s2d.rds\n" else "API: identical to api_s2d.rds\n"); quit(status = 0)
  }
  missing_ok <- setdiff(allowed, found); extra <- setdiff(found, allowed)
  if (length(missing_ok)) cat("ALLOWED CHANGE NOT PRESENT:\n", paste0("  ", missing_ok, collapse = "\n"), "\n", sep = "")
  if (length(extra)) cat("NOT ALLOWED:\n", paste0("  ", extra, collapse = "\n"), "\n", sep = "")
  cat("API: FAIL\n"); quit(status = 1)
}
stop("usage: api.R make | compare")

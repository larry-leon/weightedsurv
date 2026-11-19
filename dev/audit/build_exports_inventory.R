#!/usr/bin/env Rscript
# Builds dev/audit/exports_inventory.csv from NAMESPACE, R/, vignettes/ and man/ (read-only on all of them).
args <- commandArgs(trailingOnly = TRUE)
ns <- readLines("NAMESPACE"); ex <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
ex <- gsub("`", "", ex)
rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
rlines <- do.call(rbind, lapply(rfiles, function(f) { l <- readLines(f, warn = FALSE); data.frame(file = f, line = seq_along(l), txt = l, stringsAsFactors = FALSE) }))
code <- rlines[!grepl("^\\s*#", rlines$txt), ]
vfiles <- list.files("vignettes", pattern = "\\.Rmd$", recursive = TRUE, full.names = TRUE)
vtxt <- setNames(lapply(vfiles, readLines, warn = FALSE), vfiles)
mfiles <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
mex <- setNames(lapply(mfiles, function(f) { l <- readLines(f, warn = FALSE); i <- grep("^\\\\examples\\{", l); if (length(i)) l[i:length(l)] else character(0) }), mfiles)
esc <- function(s) gsub(".", "\\.", s, fixed = TRUE)
user_facing <- c("KM_diff","df_counting","cox_rhogamma","cumulative_rmst_bands","plotKM.band_subgroups",
                 "plot_weighted_km","KM_plot_2sample_weighted_counting","plot_weight_schemes","create_baseline_table",
                 "wt.rg.S","wlr_dhat_estimates","check_results","plot_km","get_dfcounting","km_quantile","km_quantile_table")
rows <- lapply(ex, function(e) {
  pat_def <- paste0("^`?", esc(e), "`?\\s*(<-|=)\\s*function")
  d <- rlines[grepl(pat_def, rlines$txt), ][1, ]
  pat_call <- paste0("(^|[^A-Za-z0-9._])", esc(e), "\\s*\\(")
  callers <- code[grepl(pat_call, code$txt) & !grepl(pat_def, code$txt), ]
  in_v <- names(vtxt)[vapply(vtxt, function(l) any(grepl(pat_call, l)), logical(1))]
  in_m <- names(mex)[vapply(mex, function(l) any(grepl(pat_call, l)), logical(1))]
  data.frame(export = e, file = basename(d$file), line = d$line,
             called_from_R = if (nrow(callers)) paste(unique(paste0(basename(callers$file), ":", callers$line)), collapse = "; ") else "",
             used_in_vignettes = paste(basename(in_v), collapse = "; "),
             used_in_examples = paste(sub("\\.Rd$", "", basename(in_m)), collapse = "; "),
             role = if (e %in% user_facing) "user-facing" else "helper", stringsAsFactors = FALSE)
})
inv <- do.call(rbind, rows)
if (!"--no-write" %in% args) write.csv(inv, "dev/audit/exports_inventory.csv", row.names = FALSE)
print(inv[, c("export","file","line","role")], right = FALSE)
cat("\nexports:", nrow(inv), "| user-facing:", sum(inv$role == "user-facing"), "| helper:", sum(inv$role == "helper"), "\n")
h <- inv[inv$role == "helper", ]
cat("helpers not used outside R/ (no vignette, no example use):\n"); print(h$export[h$used_in_vignettes == "" & h$used_in_examples == ""])
cat("helpers with no caller in R/:\n"); print(h$export[h$called_from_R == ""])

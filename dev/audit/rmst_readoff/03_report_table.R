#!/usr/bin/env Rscript
# 03_report_table.R -- the comparative table of REPORT_rmst_readoff_truth.md, section 3, from rmst_rows_old_new.csv
.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
x <- read.csv(file.path(.here, "rmst_rows_old_new.csv"), stringsAsFactors = FALSE)
x <- x[x$time %in% c(1, 6, 12, 24, 36, 48), ]
ord <- c(41, 42, 31, 32, 35, 37, 39, 38, 5, 10, 8); x <- x[order(match(x$cell_id, ord), x$time), ]
f <- function(v, d) ifelse(is.na(v), "—", formatC(v, format = "f", digits = d))
stored <- ifelse(x$source == "OC" & x$cell_id != 8, "R(t_j)", "R(t)")
L <- c("| k | Cell | Month | Stored CSV was against | Mean R(t) − R(t_j) | Bias vs R(t) | Bias vs R(t_j) | Cov known R(t) | Cov known R(t_j) | Cov estimated R(t) | Cov estimated R(t_j) |",
       "|---|---|---|---|---|---|---|---|---|---|---|",
       sprintf("| %d | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |", x$cell_id, x$cell, x$time, stored, f(x$gap, 5), f(x$bias_Rt, 5), f(x$bias_Rtj, 5),
               f(x$covk_Rt, 4), f(x$covk_Rtj, 4), f(x$cove_Rt, 4), f(x$cove_Rtj, 4)))
writeLines(L, file.path(.here, "report_table.md")); cat(L, sep = "\n")

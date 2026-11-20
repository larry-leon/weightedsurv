## Compares the objects of the knit of the restructured simulations article (tables and figures from
## the stored file) with those of the full-size knit of the vignette revision pass.
##   Rscript dev/audit/stage3/g6_sims_compare.R <new knit dir> <full-size knit dir>
a <- commandArgs(TRUE)
new <- readRDS(file.path(a[1], "env.rds")); old <- readRDS(file.path(a[2], "env.rds"))
for (f in list.files(a[2], "^fig-.*png$"))
  cat(sprintf("%-38s byte-identical: %s\n", f, identical(readBin(file.path(a[1], f), "raw", 1e7), readBin(file.path(a[2], f), "raw", 1e7))))
for (t in c("table_gt", "table_gt2", "table_gt3"))
  cat(sprintf("%-10s data identical: %s\n", t, identical(as.data.frame(new[[t]][["_data"]]), as.data.frame(old[[t]][["_data"]]))))
for (p in c("p_fill", "p_fill_all", "p_fill_compare"))
  cat(sprintf("%-15s plot data identical: %s\n", p, identical(as.data.frame(new[[p]]$data), as.data.frame(old[[p]]$data))))
for (o in c("n_sim", "pnull_threshold", "time_inhours", "ann_df", "ann_df_ws", "ann_df_Nonws", "min_z", "max_z"))
  cat(sprintf("%-16s identical: %s\n", o, identical(new[[o]], old[[o]])))
cat("res_demo (minimal run): rows", nrow(new$res_demo$results_sims), " number_sims", new$res_demo$number_sims, "\n")
cat("res_out rows in the new knit:", nrow(new$res_out$results_sims), " columns:", ncol(new$res_out$results_sims), "\n")

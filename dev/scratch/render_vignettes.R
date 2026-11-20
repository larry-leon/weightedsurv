#!/usr/bin/env Rscript
# render_vignettes.R -- Gate G2: rmarkdown::render() of the vignette and the articles with warnings turned
# into errors at the top level (options(warn = 2)); a file's own chunk options are left as they are.
#   Rscript dev/scratch/render_vignettes.R <lib> <outdir> [name ...]
args <- commandArgs(trailingOnly = TRUE)
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
REPO <- normalizePath(".")
lib <- normalizePath(args[1]); outdir <- args[2]
FILES <- c(examples = "vignettes/weightedsurv_examples.Rmd",
           martingale = "vignettes/articles/martingale_resampling_applications.Rmd",
           coxmeth = "vignettes/articles/weightedcox_methodology.Rmd",
           kmmeth = "vignettes/articles/weightedkm_methodology.Rmd",
           rmstmeth = "vignettes/articles/weightedrmst_methodology.Rmd",
           sims = "vignettes/articles/weightedcox_wald-vs-logrank_simulations.Rmd")
which_files <- if (length(args) > 2) args[-(1:2)] else setdiff(names(FILES), "sims")
.libPaths(c(lib, .libPaths()))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE); outdir <- normalizePath(outdir)
status <- 0
for (nm in which_files) {
  t0 <- proc.time()[["elapsed"]]
  options(warn = 2)
  res <- tryCatch({
    rmarkdown::render(file.path(REPO, FILES[[nm]]), output_dir = file.path(outdir, nm), intermediates_dir = file.path(outdir, nm),
                      envir = new.env(parent = globalenv()), quiet = TRUE)
    "OK"
  }, error = function(e) paste("FAIL:", conditionMessage(e)))
  options(warn = 0)
  cat(sprintf("%s %s (%.0f s)\n", res, FILES[[nm]], proc.time()[["elapsed"]] - t0))
  if (res != "OK") status <- 1
}
quit(status = status)

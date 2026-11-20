#!/usr/bin/env Rscript
# knit_vignettes.R -- knit the vignettes and articles with knitr::knit and keep every figure.
#   Rscript dev/scratch/knit_vignettes.R <lib> <outdir> [name ...]
# <lib>: a library holding the weightedsurv build to use. <outdir>: one sub-directory per file
# (figures fig-<chunk label>-<k>.png, the knitted .md, warnings.txt). Names: examples, martingale,
# coxmeth, kmmeth, rmstmeth, sims (default: all but sims).
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
suppressPackageStartupMessages(library(weightedsurv))
cat(sprintf("weightedsurv %s loaded from %s\n", packageVersion("weightedsurv"), dirname(getNamespaceInfo("weightedsurv", "path"))))
status <- 0
for (nm in which_files) {
  od <- file.path(outdir, nm); dir.create(od, showWarnings = FALSE)
  unlink(list.files(od, pattern = "^fig-.*\\.png$", full.names = TRUE))
  knitr::opts_chunk$set(dev = "png", dpi = 96, fig.path = file.path(od, "fig-"), error = FALSE)
  owd <- setwd(dirname(file.path(REPO, FILES[[nm]])))
  env <- new.env(parent = globalenv())
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  ok <- tryCatch({
    withCallingHandlers(
      knitr::knit(file.path(REPO, FILES[[nm]]), output = file.path(od, paste0(nm, ".md")), envir = env, quiet = TRUE),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    TRUE
  }, error = function(e) { cat(sprintf("ERROR %s: %s\n", nm, conditionMessage(e))); FALSE })
  setwd(owd)
  writeLines(warns, file.path(od, "warnings.txt"))
  saveRDS(as.list(env), file.path(od, "env.rds"))
  cat(sprintf("%s %s: %d figure file(s), %d warning(s) escaping the chunks, %.0f s\n", if (ok) "OK" else "FAIL", nm,
              length(list.files(od, pattern = "^fig-.*\\.png$")), length(warns), proc.time()[["elapsed"]] - t0))
  if (!ok) status <- 1
}
quit(status = status)

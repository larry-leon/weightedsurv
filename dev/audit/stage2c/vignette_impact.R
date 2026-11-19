#!/usr/bin/env Rscript
# =============================================================================
# vignette_impact.R -- Step V of TASK_km_fixes_stage2c.md: numbers and figures of the vignette chunks
# that pass weights, under the "before" code (START_SHA) and the "after" code (HEAD). No vignette is edited:
# each file is knitted (knitr::knit, the vignette's own code) into a scratch directory.
#
#   Rscript dev/audit/stage2c/vignette_impact.R knit <version> <lib> <outdir>
#       knits the three vignette files with weighted calls, with the package version installed in <lib>
#       first on .libPaths(); figures as PNG; after every chunk, the displayed quantities of the objects
#       that chunk created or changed are saved (snap_<file>.rds).
#   Rscript dev/audit/stage2c/vignette_impact.R compare <outdir_before> <outdir_after>
#       -> stage2c/vignette_numbers.txt (tables before -> after) and
#          stage2c/vignette_weighted_figures.pdf (each figure of an affected chunk, before | after)
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
REPO <- normalizePath(".")
S2C <- file.path(REPO, "dev", "audit", "stage2c")
FILES <- c(examples = "vignettes/weightedsurv_examples.Rmd",
           martingale = "vignettes/articles/martingale_resampling_applications.Rmd",
           coxmeth = "vignettes/articles/weightedcox_methodology.Rmd")
TT <- c(24, 60, 120)
at_t <- function(at, v, t) { i <- findInterval(t, at); ifelse(i > 0, v[pmax(i, 1)], NA_real_) }

# displayed quantities of one object (NULL when it is not one of the kinds below)
extract <- function(x) {
  if (is.list(x) && !is.data.frame(x) && !is.null(x$quantile_results) && !is.null(x$sig2_surv0)) {
    return(list(kind = "dfcount", quantiles = x$quantile_results,
                se = rbind(arm0 = at_t(x$at_points, sqrt(x$sig2_surv0), TT), arm1 = at_t(x$at_points, sqrt(x$sig2_surv1), TT)),
                surv = rbind(arm0 = at_t(x$at_points, x$surv0, TT), arm1 = at_t(x$at_points, x$surv1, TT)),
                zlogrank_text = x$zlogrank_text, z.score = x$z.score, cox_text = x$cox_results$cox_text))
  }
  if (is.list(x) && !is.null(x$fit_itt) && !is.null(x$fit_itt$dhat)) {
    f <- x$fit_itt
    return(list(kind = "band", c_alpha_band = unname(f$c_alpha_band), range = range(f$at_points),
                tab = rbind(t = TT, dhat = at_t(f$at_points, f$dhat, TT), se = at_t(f$at_points, sqrt(f$sig2_dhat), TT),
                            lower = at_t(f$at_points, f$lower, TT), upper = at_t(f$at_points, f$upper, TT),
                            sb_lower = at_t(f$at_points, f$sb_lower, TT), sb_upper = at_t(f$at_points, f$sb_upper, TT))))
  }
  if (is.list(x) && !is.null(x$rmst_text) && !is.null(x$rmst_time)) {
    n <- length(x$rmst_time)
    return(list(kind = "rmst", rmst_text = x$rmst_text, c_alpha_band = unname(x$c_alpha_band), tau = x$at_points[n],
                last = c(rmst = x$rmst_time[n], lower = x$rmst_time_lower[n], upper = x$rmst_time_upper[n],
                         sb_lower = x$rmst_time_sb_lower[n], sb_upper = x$rmst_time_sb_upper[n])))
  }
  if (is.list(x) && !is.null(x$hr_ci_asy)) {
    return(list(kind = "cox_rhogamma", z.score = x$z.score, hr_ci_asy = x$hr_ci_asy, hr_ci_star = x$hr_ci_star))
  }
  if (inherits(x, "data.table") && "analysis" %in% names(x)) return(list(kind = "table", tab = as.data.frame(x)))
  NULL
}

if (mode == "knit") {
  version <- args[2]; lib <- normalizePath(args[3]); outdir <- args[4]
  .libPaths(c(lib, .libPaths()))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE); outdir <- normalizePath(outdir)
  suppressPackageStartupMessages(library(weightedsurv))
  cat(sprintf("%s: weightedsurv loaded from %s\n", version, dirname(getNamespaceInfo("weightedsurv", "path"))))
  for (nm in names(FILES)) {
    od <- file.path(outdir, nm); dir.create(od, showWarnings = FALSE)
    SNAP <- list(); DIG <- NULL
    knitr::knit_hooks$set(snap = function(before, options, envir) {
      objs <- ls(envir)
      dg <- vapply(objs, function(o) digest::digest(get(o, envir = envir)), "")
      if (before) { DIG <<- dg; return(invisible(NULL)) }
      changed <- objs[is.na(DIG[objs]) | DIG[objs] != dg]
      ex <- Filter(Negate(is.null), lapply(setNames(changed, changed), function(o) extract(get(o, envir = envir))))
      if (length(ex)) SNAP[[options$label]] <<- ex
      invisible(NULL)
    })
    knitr::opts_chunk$set(snap = TRUE, dev = "png", dpi = 96, fig.path = file.path(od, "fig-"), error = FALSE)
    owd <- setwd(dirname(file.path(REPO, FILES[[nm]])))   # bibliography and relative paths as when the vignette is built
    env <- new.env(parent = globalenv())
    t0 <- proc.time()[["elapsed"]]
    knitr::knit(file.path(REPO, FILES[[nm]]), output = file.path(od, paste0(nm, ".md")), envir = env, quiet = TRUE)
    setwd(owd)
    saveRDS(SNAP, file.path(outdir, paste0("snap_", nm, ".rds")))
    cat(sprintf("  %s knitted in %.0f s; %d chunks with displayed quantities\n", nm, proc.time()[["elapsed"]] - t0, length(SNAP)))
  }
  quit(status = 0)
}

# chunk label -> line of its header in the Rmd (knitr numbers unnamed chunks by their position among all chunks)
chunk_lines <- function(rmd) {
  L <- readLines(file.path(REPO, rmd)); h <- grep("^```\\{r", L)
  lab <- vapply(seq_along(h), function(i) {
    m <- regmatches(L[h[i]], regexec("^```\\{r[ ,]*([A-Za-z0-9_.-]+)?\\s*[,}]", L[h[i]]))[[1]]
    if (length(m) >= 2 && nzchar(m[2]) && !grepl("=", m[2])) m[2] else paste0("unnamed-chunk-", i)
  }, "")
  setNames(h, lab)
}

if (mode == "compare") {
  B <- args[2]; A <- args[3]
  sink(file.path(S2C, "vignette_numbers.txt"), split = TRUE)
  cat("Vignette numbers, before (START_SHA) -> after (HEAD). Times for SE / band read-off:", paste(TT, collapse = ", "), "months.\n")
  fig_list <- list()
  fmt <- function(x) if (is.numeric(x)) paste(formatC(x, digits = 4, format = "g"), collapse = " / ") else paste(x, collapse = " / ")
  for (nm in names(FILES)) {
    sb <- readRDS(file.path(B, paste0("snap_", nm, ".rds"))); sa <- readRDS(file.path(A, paste0("snap_", nm, ".rds")))
    cl <- chunk_lines(FILES[[nm]])
    cat(sprintf("\n==== %s (%s) ====\n", nm, FILES[[nm]]))
    for (lab in union(names(sb), names(sa))) {
      for (o in union(names(sb[[lab]]), names(sa[[lab]]))) {
        b <- sb[[lab]][[o]]; a <- sa[[lab]][[o]]
        same <- identical(b, a)
        cat(sprintf("\n-- chunk '%s' (line %s), object %s [%s]: %s\n", lab, cl[lab], o, (if (is.null(a)) b else a)$kind,
                    if (same) "IDENTICAL before/after" else "CHANGED"))
        if (same) next
        for (f in setdiff(names(a), "kind")) {
          if (identical(b[[f]], a[[f]])) { cat(sprintf("   %-14s unchanged: %s\n", f, if (is.data.frame(a[[f]]) || is.matrix(a[[f]])) "(table)" else fmt(a[[f]]))); next }
          if (is.data.frame(a[[f]]) || is.matrix(a[[f]])) {
            cat(sprintf("   %s before:\n", f)); print(b[[f]], digits = 4); cat(sprintf("   %s after:\n", f)); print(a[[f]], digits = 4)
          } else cat(sprintf("   %-14s %s  ->  %s\n", f, fmt(b[[f]]), fmt(a[[f]])))
        }
      }
    }
    # figures of every chunk that passes weights (figure files present under both versions)
    figs_b <- list.files(file.path(B, nm), pattern = "^fig-.*\\.png$"); figs_a <- list.files(file.path(A, nm), pattern = "^fig-.*\\.png$")
    for (fg in intersect(figs_b, figs_a)) {
      lab <- sub("^fig-(.*)-[0-9]+\\.png$", "\\1", fg)
      fig_list[[length(fig_list) + 1]] <- list(file = nm, label = lab, line = cl[lab], fig = fg,
        b = file.path(B, nm, fg), a = file.path(A, nm, fg),
        same = identical(unname(tools::md5sum(file.path(B, nm, fg))), unname(tools::md5sum(file.path(A, nm, fg)))))
    }
  }
  cat("\n==== figures (PNG files byte-identical before/after?) ====\n")
  for (f in fig_list) cat(sprintf("  %s line %s %s: %s\n", f$file, f$line, f$fig, if (f$same) "identical" else "DIFFERS"))
  sink()
  # side-by-side sheet of the figures that differ
  grDevices::pdf(file.path(S2C, "vignette_weighted_figures.pdf"), width = 14, height = 6)
  for (f in Filter(function(f) !f$same, fig_list)) {
    graphics::par(mar = c(0, 0, 2.5, 0), mfrow = c(1, 2), oma = c(0, 0, 2, 0))
    for (side in c("b", "a")) {
      im <- png::readPNG(f[[side]]); asp <- dim(im)[1] / dim(im)[2]
      graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, asp), asp = 1, xaxs = "i", yaxs = "i")
      graphics::rasterImage(im, 0, 0, 1, asp)
      graphics::title(main = if (side == "b") "before (START_SHA, Stage 2b code)" else "after (HEAD, Stage 2c code)", cex.main = 1)
    }
    graphics::mtext(sprintf("%s, chunk at line %s (%s)", FILES[[f$file]], f$line, f$fig), outer = TRUE, cex = 1.1, font = 2)
  }
  grDevices::dev.off()
  cat(sprintf("vignette_weighted_figures.pdf: %d pages (figures that differ)\n", sum(!vapply(fig_list, `[[`, TRUE, "same"))))
  quit(status = 0)
}
stop("usage: vignette_impact.R knit <version> <lib> <outdir> | compare <before> <after>")

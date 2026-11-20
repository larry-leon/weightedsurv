#!/usr/bin/env Rscript
# =============================================================================
# vignette_impact.R -- Step V of TASK_km_fixes_stage2d.md: numbers and figures of the vignette chunks
# that pass weights, under the "before" code (START_SHA, Stage 2c) and the "after" code (HEAD, Stage 2d).
# No vignette is edited: each file is knitted (knitr::knit, the vignette's own code) into a scratch directory.
# Copy of dev/audit/stage2c/vignette_impact.R, changed as follows: output goes to stage2d/
# (vignette_numbers_2d.txt, vignette_figures_2d.pdf); the displayed quantities recorded for a df_counting()
# result include the log-rank score, its variance and z; those of a cox_rhogamma() result include both
# standard errors and both intervals; a wlr_cumulative() result is recorded; and, after each file is knitted,
# a supplement re-runs, in the knitted environment, the cox_rhogamma() calls of the chunks that overwrite
# their result (one object per scheme), so that every displayed scheme is recorded.
#
#   Rscript dev/audit/stage2d/vignette_impact.R knit <version> <lib> <outdir>
#   Rscript dev/audit/stage2d/vignette_impact.R compare <outdir_before> <outdir_after>
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
REPO <- normalizePath(".")
S2C <- file.path(REPO, "dev", "audit", "stage2d")
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
                zlogrank_text = x$zlogrank_text, z.score = x$z.score, lr = x$lr, sig2_lr = x$sig2_lr,
                z.score_stratified = x$z.score_stratified, cox_text = x$cox_results$cox_text,
                logrank_text_survdiff = x$logrank_results$logrank_text))
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
    return(list(kind = "cox_rhogamma", z.score = x$z.score, zlogrank_text = x$zlogrank_text, bhat = x$fit$bhat,
                sig_bhat_asy = x$fit$sig_bhat_asy, sig_bhat_star = x$fit$sig_bhat_star, bhat_debiased = x$fit$bhat_debiased,
                hr_ci_asy = x$hr_ci_asy, hr_ci_star = x$hr_ci_star, cox_text_asy = x$cox_text_asy, cox_text_star = x$cox_text_star))
  }
  if (is.list(x) && !is.data.frame(x) && identical(names(x), c("z.score", "time", "score", "sig2.score")) && length(x$time) > 1) {
    n <- length(x$time)
    return(list(kind = "wlr_cumulative", last = c(score = x$score[n], sig2.score = x$sig2.score[n], z.score = x$z.score[n]),
                z_range = range(x$z.score, na.rm = TRUE), first_time_z_above_1.96 = x$time[which(x$z.score > stats::qnorm(0.975))[1]],
                time_of_max_z = x$time[which.max(x$z.score)], n_grid_points_z_above_5 = sum(x$z.score > 5, na.rm = TRUE),
                z_range_from_80 = range(x$z.score[x$time >= 80], na.rm = TRUE),
                tab = rbind(t = TT, score = at_t(x$time, x$score, TT), sig2.score = at_t(x$time, x$sig2.score, TT), z.score = at_t(x$time, x$z.score, TT))))
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
    # supplement: the cox_rhogamma() calls whose result the chunk overwrites, one object per scheme
    sup <- list()
    cxs <- function(dc, scheme, sp = list(rho = 0, gamma = 0.5)) cox_rhogamma(dfcount = dc, scheme = scheme, scheme_params = sp, draws = 1000, verbose = FALSE)
    if (nm == "examples" && exists("dfcount_rotterdam_wtd", envir = env)) {
      dc <- get("dfcount_rotterdam_wtd", envir = env)
      # the chunk at line 711 re-creates dfcount_rotterdam_wtd with draws = 5000; the log-rank and Cox inputs are the same
      sup[["cox FH(0,0)"]] <- extract(cxs(dc, "fh", list(rho = 0, gamma = 0)))
      sup[["cox FH(0,1)"]] <- extract(cxs(dc, "fh", list(rho = 0, gamma = 1)))
      sup[["cox FH_exp2"]] <- extract(cox_rhogamma(dfcount = dc, scheme = "fh_exp2", draws = 1000, verbose = FALSE))
      sup[["cox MB(12)"]] <- extract(cxs(dc, "MB", list(mb_tstar = 12)))
      sup[["cox zero_one(6)"]] <- extract(cxs(dc, "custom_time", list(t.tau = 6, w0.tau = 0, w1.tau = 1)))
      SNAP[["supplement-624 (cox_rhogamma on dfcount_rotterdam_wtd, draws = 1000, each scheme)"]] <- sup
    }
    if (nm == "examples" && exists("df_CTC", envir = env) && exists("dfcount_CTC_wtd", envir = env)) {
      dctc <- get("df_CTC", envir = env); dcc <- get("dfcount_CTC_wtd", envir = env)
      cw <- wlr_cumulative(dcc, scheme_params = list(rho = 0, gamma = 0), scheme = "fh", return_cumulative = TRUE)
      tp <- cw$time[which.max(cw$z.score)]
      w <- dctc$sw.weights; ev <- dcc$delta == 1 & dcc$time <= tp
      SNAP[["supplement-502 (cross-trial data up to the time of the largest cumulative z)"]] <- list(ctc = list(kind = "note",
        time_of_max_z = tp, n_by_arm_0_1 = as.numeric(table(dcc$z)), events_to_that_time_arm_0_1 = c(sum(ev & dcc$z == 0), sum(ev & dcc$z == 1)),
        weighted_events_to_that_time_arm_0_1 = c(sum(dcc$w_hat[ev & dcc$z == 0]), sum(dcc$w_hat[ev & dcc$z == 1])),
        weight_range = range(w), weight_range_arm1 = range(w[dctc$treat == 1]), effective_n_arm_0_1 = c(sum(w[dctc$treat == 0])^2 / sum(w[dctc$treat == 0]^2), sum(w[dctc$treat == 1])^2 / sum(w[dctc$treat == 1]^2))))
    }
    if (nm == "examples" && exists("df_rotterdam", envir = env)) {
      dr <- get("df_rotterdam", envir = env); w <- dr$sw.weights
      SNAP[["supplement-299 (Rotterdam weights)"]] <- list(rot = list(kind = "note", n_by_arm_0_1 = as.numeric(table(dr$treat)), weight_range = range(w),
        weight_range_arm1 = range(w[dr$treat == 1]), effective_n_arm_0_1 = c(sum(w[dr$treat == 0])^2 / sum(w[dr$treat == 0]^2), sum(w[dr$treat == 1])^2 / sum(w[dr$treat == 1]^2))))
    }
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
  sink(file.path(S2C, "vignette_numbers_2d.txt"), split = TRUE)
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
  grDevices::pdf(file.path(S2C, "vignette_figures_2d.pdf"), width = 14, height = 6)
  for (f in Filter(function(f) !f$same, fig_list)) {
    graphics::par(mar = c(0, 0, 2.5, 0), mfrow = c(1, 2), oma = c(0, 0, 2, 0))
    for (side in c("b", "a")) {
      im <- png::readPNG(f[[side]]); asp <- dim(im)[1] / dim(im)[2]
      graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, asp), asp = 1, xaxs = "i", yaxs = "i")
      graphics::rasterImage(im, 0, 0, 1, asp)
      graphics::title(main = if (side == "b") "before (START_SHA, Stage 2c code)" else "after (HEAD, Stage 2d code)", cex.main = 1)
    }
    graphics::mtext(sprintf("%s, chunk at line %s (%s)", FILES[[f$file]], f$line, f$fig), outer = TRUE, cex = 1.1, font = 2)
  }
  grDevices::dev.off()
  cat(sprintf("vignette_figures_2d.pdf: %d pages (figures that differ)\n", sum(!vapply(fig_list, `[[`, TRUE, "same"))))
  quit(status = 0)
}
stop("usage: vignette_impact.R knit <version> <lib> <outdir> | compare <before> <after>")

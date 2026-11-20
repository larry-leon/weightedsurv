#!/usr/bin/env Rscript
# gates_appendix.R -- Gates G3 (time, code grep) and G4 of TASK_estimated_weights_appendix.md.
# The chunks of the new article are run with knitr::knit() in a fresh environment; the by-hand and
# package values of its section 8 are then read from that environment.
#   Rscript dev/scratch/gates_appendix.R <lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(normalizePath(args[1]), .libPaths()))
f <- normalizePath("vignettes/articles/estimated_weights_methodology.Rmd")
e <- new.env(parent = globalenv())
td <- tempfile("ew"); dir.create(td); owd <- setwd(td)
t0 <- proc.time()[["elapsed"]]
invisible(knitr::knit(f, output = file.path(td, "out.md"), envir = e, quiet = TRUE))
el <- proc.time()[["elapsed"]] - t0
setwd(owd)
cat(sprintf("%s G3 (time): chunks of the article run in %.1f s (limit 60)\n", if (el < 60) "PASS" else "FAIL", el))
# code lines only (inside fenced chunks)
ln <- readLines(f); fence <- cumsum(grepl("^```", ln)) %% 2 == 1
code <- ln[fence & !grepl("^```", ln)]
hits <- grep("replicate\\(|boot|for \\(b in|mclapply|future|glm\\(", code, value = TRUE)
cat(sprintf("%s G3 (code): %d code lines; lines matching replicate( | boot | for (b in | mclapply | future | glm( : %d; calls to ps_weights(): %d\n",
            if (length(hits) == 0) "PASS" else "FAIL", length(code), length(hits), length(grep("ps_weights\\(", code))))
cmp <- e$cmp
cat("G4: by hand / package / difference\n")
for (i in seq_len(nrow(cmp))) cat(sprintf("   %-45s %.15f  %.15f  %.2e\n", cmp$quantity[i], cmp$by_hand[i], cmp$package[i], cmp$difference[i]))
ok <- all(abs(cmp$difference) < 1e-8)
cat(sprintf("%s G4: largest absolute difference %.2e (limit 1e-8)\n", if (ok) "PASS" else "FAIL", max(abs(cmp$difference))))
cat(sprintf("   also in the article: n = %d, treated %d; ESS ratio treated %.4f, control %.4f; contributions by formula against cox_rhogamma_resample() with unit multipliers: %.1e\n",
            nrow(e$df_rott), sum(e$df_rott$treat == 1), e$dc_est$ess$ratio[1], e$dc_est$ess$ratio[2],
            max(abs(-e$unit$score_star - e$sc_bhat$psi), abs(-e$unit$score_star_null - e$sc_zero$psi))))
cat(sprintf("   estimates identical under the two options: dhat %s, bhat %s\n", identical(e$fit_est$dhat, e$fit_known$dhat), identical(e$cox_est$fit$bhat, e$cox_known$fit$bhat)))
print(e$dec)
# Lemmas 1 and 2 against central differences in one subject's weight (five subjects of each arm), from exported functions
with(as.list(e), {
  h <- 1e-6; set.seed(1); ks <- c(sample(which(A == 1), 5), sample(which(A == 0), 5))
  dfs <- data.frame(tte = U, event = D, treat = A)
  S_arm <- function(ww, a) { f <- KM_diff(transform(dfs, w = ww), "tte", "event", "treat", weight.name = "w", at_points = t_star); if (a == 1) f$surv1 else f$surv0 }
  r1 <- sapply(ks, function(k) { up <- w; up[k] <- up[k] + h; dn <- w; dn[k] <- dn[k] - h; a <- A[k]; kk <- if (a == 1) k1 else k0
    c((S_arm(up, a) - S_arm(dn, a)) / (2 * h), -kk$S * kk$psi[k]) })
  r2 <- sapply(ks, function(k) { up <- w; up[k] <- up[k] + h; dn <- w; dn[k] <- dn[k] - h
    c((cox_score_rhogamma(bhat, U, D, A, up) - cox_score_rhogamma(bhat, U, D, A, dn)) / (2 * h), -sc_bhat$psi[k] / w[k]) })
  cat(sprintf("   Lemma 1 (dS^a/dw_i = -S^a psi_i), 10 subjects: largest relative difference from the central difference %.1e\n", max(abs(r1[1, ] - r1[2, ]) / abs(r1[2, ]))))
  cat(sprintf("   Lemma 2 (dU/dw_j = -psi_j / w_j at bhat), 10 subjects: largest relative difference from the central difference %.1e\n", max(abs(r2[1, ] - r2[2, ]) / abs(r2[2, ]))))
  cat(sprintf("   I_w(bhat) by the formula of section 4.4 against cox_score_rhogamma()$i_bhat: %.1e\n",
              abs(sc_bhat$info - cox_score_rhogamma(bhat, U, D, A, w, score_only = FALSE)$i_bhat)))
})
quit(status = if (ok && el < 60 && length(hits) == 0) 0 else 1)

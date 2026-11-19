# band_draws20_check.R -- plotKM.band_subgroups() defaults to draws = 20, so with weights its pointwise SE
# (and the band standardisation) is a 20-draw resampling variance. Compares it with the analytic (draws = 0)
# SE on the vignette's Rotterdam propensity-weighted data, at HEAD. Writes band_draws20_check.txt.
suppressMessages(pkgload::load_all(".", quiet = TRUE)); grDevices::pdf(NULL)
r <- within(survival::rotterdam, { rfstime <- ifelse(recur == 1, rtime, dtime); status <- pmax(recur, death)
  ignore <- (recur == 0 & death == 1 & rtime < dtime); status2 <- ifelse(recur == 1 | ignore, recur, death)
  grade3 <- ifelse(grade == "3", 1, 0); treat <- hormon; event <- status2; tte <- rfstime / 30.4375 })
d <- subset(r, nodes > 0)
ps <- glm(treat ~ age + meno + size + grade3 + nodes + pgr + chemo + er, data = d, family = "binomial")$fitted
p0 <- glm(treat ~ 1, family = "binomial", data = d)$fitted
d$sw.weights <- ifelse(d$treat == 1, p0 / ps, (1 - p0) / (1 - ps))
b <- plotKM.band_subgroups(df = d, tte.name = "tte", treat.name = "treat", event.name = "event", weight.name = "sw.weights",
                           tau_add = 42, by.risk = 12, draws.band = 1000, qtau = 0.025, show_resamples = FALSE)$fit_itt
k0 <- KM_diff(d, "tte", "event", "treat", weight.name = "sw.weights", at_points = b$at_points)
i <- findInterval(c(12, 24, 36, 60, 90, 120), b$at_points)
out <- data.frame(t = c(12, 24, 36, 60, 90, 120), se_draws20 = sqrt(b$sig2_dhat[i]), se_analytic = sqrt(k0$sig2_dhat[i]))
out$ratio <- out$se_draws20 / out$se_analytic
r_all <- sqrt(b$sig2_dhat / k0$sig2_dhat)[k0$sig2_dhat > 0]
sink("dev/audit/stage2c/band_draws20_check.txt", split = TRUE)
print(out, digits = 4); cat(sprintf("over the band grid: ratio quartiles %s, range %.3f to %.3f\n",
  paste(sprintf("%.3f", quantile(r_all, c(.25, .5, .75))), collapse = " / "), min(r_all), max(r_all))); sink()

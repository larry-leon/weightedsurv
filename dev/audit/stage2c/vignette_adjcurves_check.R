# vignette_adjcurves_check.R <lib> -- the SE comparison of weightedsurv_examples.Rmd (lines 702-768) under the
# package version installed in <lib>: package SE / adjustedCurves (iptw_km, conf_int) SE and package SE /
# survfit(weights, robust = TRUE) SE, at the event times, for draws = 0 and draws = 5000 (the vignette's df_mine).
lib <- commandArgs(TRUE)[1]; .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(survival); library(weightedsurv); library(adjustedCurves) }); grDevices::pdf(NULL)
g <- within(rotterdam, { rfstime <- ifelse(recur == 1, rtime, dtime); time_months <- rfstime / 30.4375; status <- pmax(recur, death)
  ignore <- (recur == 0 & death == 1 & rtime < dtime); status2 <- ifelse(recur == 1 | ignore, recur, death)
  grade3 <- ifelse(grade == "3", 1, 0); treat <- hormon; event <- status2; tte <- time_months })
d <- subset(g, nodes > 0)
ps <- glm(treat ~ age + meno + size + grade3 + nodes + pgr + chemo + er, data = d, family = "binomial")$fitted
p0 <- glm(treat ~ 1, family = "binomial", data = d)$fitted
d$sw.weights <- ifelse(d$treat == 1, p0 / ps, (1 - p0) / (1 - ps))
d$hormon2 <- factor(d$hormon, labels = c("No", "Yes"))
pm <- glm(hormon2 ~ age + meno + size + grade + nodes + pgr + chemo + er, data = d, family = "binomial")
iptw <- suppressWarnings(adjustedsurv(data = d, variable = "hormon2", ev_time = "time_months", event = "status", method = "iptw_km", treatment_model = pm, conf_int = TRUE))
cat("weightedsurv from", dirname(getNamespaceInfo("weightedsurv", "path")), "\n")
for (dr in c(0, 5000)) {
  m <- suppressWarnings(get_dfcounting(df = d, tte.name = "tte", event.name = "event", treat.name = "treat", arms = c("treat", "control"), by.risk = 24, weight.name = "sw.weights", draws = dr))
  for (a in 0:1) {
    ch <- subset(iptw$adj, group == c("No", "Yes")[a + 1]); ch <- ch[ch$se > 0 & is.finite(ch$se), ]
    idv <- m[[paste0("idv", a, ".check")]]; tt <- m$at_points[idv]; se <- sqrt(m[[paste0("sig2_surv", a)]])[idv]
    ref <- ch$se[findInterval(tt, ch$time)]
    x <- d[d$treat == a, ]; rb <- summary(survfit(aeqSurv(Surv(x$tte, x$event)) ~ 1, weights = x$sw.weights, robust = TRUE), times = tt, extend = TRUE)$std.err
    k <- is.finite(ref) & ref > 0 & rb > 0
    cat(sprintf("draws = %4d, arm %d: pkg/adjustedCurves SE median %.3f (range %.3f-%.3f); pkg/survfit robust median %.4f (range %.4f-%.4f)\n",
                dr, a, median(se[k] / ref[k]), min(se[k] / ref[k]), max(se[k] / ref[k]), median(se[k] / rb[k]), min(se[k] / rb[k]), max(se[k] / rb[k])))
  }
}

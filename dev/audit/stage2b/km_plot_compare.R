#!/usr/bin/env Rscript
# =============================================================================
# km_plot_compare.R <before|after> -- Kaplan-Meier plot comparison against survfit() on GBSG
# (TASK_km_fixes_stage2.md, Step P)
#
# Writes dev/audit/stage2/km_plots_<mode>.pdf  (KM-curve sheets, then difference-curve sheets)
#        dev/audit/stage2/km_plot_compare_<mode>.csv
#
# "Drawn" = what the package's plotting code puts on the device:
#   before, or after without the B10 helper: at_points[idv], surv[idv] joined with type = "s";
#            the band is polygon() over the same vertices, so it is linear between them.
#   after with the B10 helper (internal km_step_vertices()): its step vertices for the curve and
#            the band (both step functions).
# Differences are evaluated over [0, last time] at every survfit jump and at the midpoints
# between consecutive survfit times, wherever the package draws something (first_drawn_x
# shows where drawing starts). Band comparisons skip points where survfit's limit is NA (S = 0).
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else stop("usage: km_plot_compare.R before|after")
stopifnot(mode %in% c("before", "after"))
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
S2 <- file.path("dev", "audit", "stage2")
ns <- asNamespace("weightedsurv")
HELPER <- mode == "after" && exists("km_step_vertices", envir = ns, inherits = FALSE)
cat(sprintf("mode %s; B10 vertex helper %s\n", mode, if (HELPER) "used" else "not used"))

g <- survival::gbsg
gb <- data.frame(tte = g$rfstime / 30.4375, event = g$status, treat = g$hormon,
                 er = g$er, pgr = g$pgr, meno = g$meno, grade = g$grade, nodes = g$nodes, size = g$size, age = g$age)
SG <- c("er <= 0", "er > 0", "er <= 5", "er > 5", "er <= 20", "er > 20", "er <= 50", "er > 50", "er > 100",
        "er <= 0 & meno == 0", "er <= 0 & meno == 1", "er <= 0 & grade == 3", "er <= 0 & nodes > 3",
        "er <= 0 & pgr <= 0", "er <= 0 & age <= 50", "er > 0 & pgr > 0", "er > 0 & grade == 3", "er > 0 & age <= 40",
        "pgr <= 0", "pgr > 0", "pgr <= 20", "pgr > 20",
        "meno == 0", "meno == 1", "grade == 1", "grade == 2", "grade == 3",
        "nodes <= 3", "nodes > 3", "nodes > 9", "size <= 20", "size > 20", "size > 50",
        "age <= 40", "age > 40 & age <= 60", "age > 60", "(all patients)")
sg_expr <- function(s) if (s == "(all patients)") "TRUE" else s
COL <- c("0" = "black", "1" = "blue"); LIGHT <- c("0" = "grey75", "1" = "lightblue")

capture <- function(expr) {
  msgs <- character(0); ws <- character(0)
  val <- tryCatch(withCallingHandlers(expr,
           message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
           warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }),
         error = function(e) structure(conditionMessage(e), class = "cap_error"))
  list(value = if (inherits(val, "cap_error")) NULL else val, error = if (inherits(val, "cap_error")) unclass(val) else NULL,
       messages = msgs, warnings = ws)
}

# what the package draws for one arm
drawn <- function(dc, arm) {
  at <- dc$at_points; S <- if (arm == 0) dc$surv0 else dc$surv1
  se <- sqrt(if (arm == 0) dc$sig2_surv0 else dc$sig2_surv1); idv <- if (arm == 0) dc$idv0 else dc$idv1
  if (HELPER) {
    v <- get("km_step_vertices", envir = ns)(at_points = at, S = S, idv = idv, se = se, conf_level = 0.95, time.zero = 0)
    return(list(x = v$x, y = v$y, bx = v$x, lower = v$lower, upper = v$upper, band = "step"))
  }
  z <- stats::qnorm(0.975)
  x <- at[idv]; y <- S[idv]; s <- se[idv]
  list(x = x, y = y, bx = x, lower = exp(log(y) - z * s / y), upper = exp(log(y) + z * s / y), band = "linear")
}
step_eval <- function(x, y, t) { i <- findInterval(t, x); out <- rep(NA_real_, length(t)); out[i > 0] <- y[i[i > 0]]; out }
lin_eval <- function(x, y, t) {
  ok <- is.finite(y); out <- rep(NA_real_, length(t))
  if (sum(ok) >= 2) { r <- t >= min(x) & t <= max(x); out[r] <- stats::approx(x[ok], y[ok], t[r], ties = "ordered")$y }
  out
}
mx <- function(v) { v <- v[is.finite(v)]; if (length(v)) max(v) else NA_real_ }

rows <- list()
pdf_file <- file.path(S2, sprintf("km_plots_%s.pdf", mode))
grDevices::pdf(pdf_file, width = 11, height = 8.5)
for (s in SG) {
  sub <- gb[with(gb, eval(parse(text = sg_expr(s)))), ]
  ne <- tapply(sub$event, factor(sub$treat, 0:1), sum); nn <- table(factor(sub$treat, 0:1))
  title_txt <- sprintf("%s  |  control n = %d, events %d  |  treated n = %d, events %d", s, nn[1], ne[1], nn[2], ne[2])
  dcr <- capture(df_counting(sub, "tte", "event", "treat", by.risk = 12))
  dc <- dcr$value
  graphics::par(mfrow = c(2, 2), oma = c(0, 0, 2, 0), mar = c(4, 4, 2, 1))
  # panel 1: package
  if (is.null(dc)) {
    graphics::plot.new(); graphics::title(main = "plot_weighted_km()"); graphics::text(0.5, 0.5, paste("df_counting error:", dcr$error), cex = 0.7)
  } else {
    pr <- capture(plot_weighted_km(dc, conf.int = TRUE, show.cox = FALSE, show.med = FALSE))
    if (length(pr$messages) || !is.null(pr$error)) {
      if (grDevices::dev.cur() > 1 && !length(grDevices::recordPlot()[[1]])) graphics::plot.new()
      graphics::mtext(paste("package:", paste(c(pr$messages, pr$error), collapse = " ")), side = 3, line = 0, cex = 0.6, col = "red")
    }
    graphics::title(main = "plot_weighted_km(conf.int = TRUE)", cex.main = 0.9)
  }
  # panel 2: survfit
  fit <- survival::survfit(survival::Surv(tte, event) ~ treat, data = sub, conf.type = "log")
  graphics::plot(fit, conf.int = TRUE, mark.time = TRUE, col = COL, ylim = c(-0.125, 1), xlim = c(0, max(sub$tte)),
                 xlab = "time", ylab = "proportion surviving", main = "survfit(conf.int = TRUE)", cex.main = 0.9)
  # panels 3 and 4: overlay
  for (xl in list(c(0, max(sub$tte)), c(0, 12))) {
    graphics::plot(NA, xlim = xl, ylim = c(0, 1), xlab = "time", ylab = "S(t)",
                   main = if (xl[2] == 12) "overlay, left tail (0-12)" else "overlay: survfit thick light, package vertices thin dark", cex.main = 0.9)
    for (a in 0:1) {
      fa <- survival::survfit(survival::Surv(tte, event) ~ 1, data = sub[sub$treat == a, ])
      graphics::lines(c(0, fa$time), c(1, fa$surv), type = "s", col = LIGHT[as.character(a)], lwd = 7)
      if (!is.null(dc)) { v <- drawn(dc, a); graphics::lines(v$x, v$y, type = "s", col = COL[as.character(a)], lwd = 1) }
    }
  }
  graphics::mtext(title_txt, outer = TRUE, cex = 0.9, font = 2)

  # numbers
  for (a in 0:1) {
    d1 <- sub[sub$treat == a, ]
    fe <- if (any(d1$event == 1)) min(d1$tte[d1$event == 1]) else NA_real_
    row <- data.frame(subgroup = s, arm = a, n = nrow(d1), events = sum(d1$event), first_event_time = fe,
                      censored_before_first_event = if (is.na(fe)) sum(d1$event == 0) else sum(d1$event == 0 & d1$tte < fe),
                      first_drawn_x = NA_real_, first_drawn_y = NA_real_, band_upper_max = NA_real_,
                      max_abs_vertex_diff_curve = NA_real_, max_abs_vertex_diff_band = NA_real_, max_abs_dhat_diff = NA_real_,
                      status = if (is.null(dc)) paste("df_counting error:", dcr$error) else "ok", stringsAsFactors = FALSE)
    if (!is.null(dc) && nrow(d1) > 0) {
      v <- drawn(dc, a)
      fa <- survival::survfit(survival::Surv(tte, event) ~ 1, data = d1, conf.type = "log")
      sx <- c(0, fa$time); sy <- c(1, fa$surv); sl <- c(1, fa$lower); su <- c(1, fa$upper)
      last <- max(d1$tte)
      tt <- sort(unique(c(sx, (utils::head(sx, -1) + utils::tail(sx, -1)) / 2))); tt <- tt[tt <= last]
      ref_S <- step_eval(sx, sy, tt); got_S <- step_eval(v$x, v$y, tt)
      ref_l <- step_eval(sx, sl, tt); ref_u <- step_eval(sx, su, tt)
      if (v$band == "step") { got_l <- step_eval(v$bx, v$lower, tt); got_u <- step_eval(v$bx, v$upper, tt) }
      else { got_l <- lin_eval(v$bx, v$lower, tt); got_u <- lin_eval(v$bx, v$upper, tt) }
      row$first_drawn_x <- v$x[1]; row$first_drawn_y <- v$y[1]
      row$band_upper_max <- mx(v$upper)
      row$max_abs_vertex_diff_curve <- mx(abs(got_S - ref_S))
      row$max_abs_vertex_diff_band <- mx(c(abs(got_l - ref_l), abs(got_u - ref_u)))
    }
    rows[[length(rows) + 1]] <- row
  }
}

# difference-curve sheets: plotKM.band_subgroups, four subgroups at a time, beside survfit
SGC <- c("red", "darkgreen", "purple", "orange")
dh <- list()
for (k in seq(1, length(SG), by = 4)) {
  labs <- SG[k:min(k + 3, length(SG))]
  graphics::par(mfrow = c(1, 2), oma = c(0, 0, 2, 0), mar = c(5, 5, 2, 1))
  br <- capture(plotKM.band_subgroups(gb, "tte", "event", "treat", sg_labels = vapply(labs, sg_expr, ""),
                                      sg_colors = SGC[seq_along(labs)], draws.band = 0))
  if (is.null(br$value)) {
    graphics::plot.new(); graphics::text(0.5, 0.5, paste("plotKM.band_subgroups error:", br$error), cex = 0.7)
    graphics::plot.new()
    for (s in labs) dh[[s]] <- NA_real_
  } else {
    xp <- br$value$xpoints
    ref <- sapply(labs, function(s) {
      sub <- gb[with(gb, eval(parse(text = sg_expr(s)))), ]
      sv <- sapply(0:1, function(a) summary(survival::survfit(survival::Surv(tte, event) ~ 1, data = sub[sub$treat == a, ]),
                                             times = xp, extend = TRUE)$surv)
      sv[, 2] - sv[, 1]
    })
    ref <- matrix(ref, ncol = length(labs))
    for (i in seq_along(labs)) dh[[labs[i]]] <- mx(abs(br$value$Dhat_subgroups[, i] - ref[, i]))
    graphics::matplot(xp, ref, type = "s", lty = 1, lwd = 2, col = SGC[seq_along(labs)], ylim = graphics::par("usr")[3:4],
                      xlab = "Months", ylab = "S1 - S0 (survfit)", main = "survfit differences at xpoints", cex.main = 0.9)
    graphics::abline(h = 0, lty = 2)
    graphics::legend("bottomleft", legend = labs, col = SGC[seq_along(labs)], lty = 1, lwd = 2, bty = "n", cex = 0.7)
  }
  graphics::mtext(sprintf("plotKM.band_subgroups(draws.band = 0): %s", paste(labs, collapse = " | ")), outer = TRUE, cex = 0.8, font = 2)
}
grDevices::dev.off()

res <- do.call(rbind, rows)
res$max_abs_dhat_diff <- unlist(dh[res$subgroup])
num <- c("first_event_time", "first_drawn_x", "first_drawn_y", "band_upper_max", "max_abs_vertex_diff_curve", "max_abs_vertex_diff_band", "max_abs_dhat_diff")
res[num] <- lapply(res[num], signif, 6)
utils::write.csv(res, file.path(S2, sprintf("km_plot_compare_%s.csv", mode)), row.names = FALSE)
cat(sprintf("wrote %s (%d subgroups) and km_plot_compare_%s.csv (%d rows)\n", pdf_file, length(SG), mode, nrow(res)))

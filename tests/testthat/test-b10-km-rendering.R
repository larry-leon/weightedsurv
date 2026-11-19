# B10: the KM curve and its band are drawn as step functions from time zero, equal to survfit's.

step_at <- function(x, y, t) { i <- findInterval(t, x); y[i] }
gbsg_er0 <- function() {
  g <- survival::gbsg[survival::gbsg$er <= 0, ]
  data.frame(tte = g$rfstime / 30.4375, event = g$status, treat = g$hormon)
}

test_that("vertices: step function from (0, 1) equal to survfit's curve and log limits", {
  d <- gbsg_er0()
  dc <- df_counting(d, "tte", "event", "treat", by.risk = 12)
  for (a in 0:1) {
    S <- if (a == 0) dc$surv0 else dc$surv1
    se <- sqrt(if (a == 0) dc$sig2_surv0 else dc$sig2_surv1)
    v <- km_step_vertices(dc$at_points, S, if (a == 0) dc$idv0 else dc$idv1, se = se, conf_level = 0.95)
    expect_equal(c(v$x[1], v$y[1]), c(0, 1))
    fa <- survival::survfit(survival::Surv(tte, event) ~ 1, data = d[d$treat == a, ], conf.type = "log")
    sx <- c(0, fa$time)
    tt <- sort(unique(c(sx, (utils::head(sx, -1) + utils::tail(sx, -1)) / 2)))
    tt <- tt[tt <= max(d$tte[d$treat == a])]
    expect_equal(step_at(v$x, v$y, tt), step_at(sx, c(1, fa$surv), tt), tolerance = 1e-10)
    expect_equal(step_at(v$x, v$lower, tt), step_at(sx, c(1, fa$lower), tt), tolerance = 1e-10)
    expect_equal(step_at(v$x, v$upper, tt), step_at(sx, c(1, fa$upper), tt), tolerance = 1e-10)
    expect_lte(max(v$upper), 1)
    expect_equal(utils::tail(v$x, 1), max(d$tte[d$treat == a]))
  }
})

test_that("weighted curve vertices equal survfit(weights =)", {
  d <- ws_sim(150, 61, weighted = TRUE)
  dc <- df_counting(d, "time", "status", "treat", weight.name = "w")
  v <- km_step_vertices(dc$at_points, dc$surv1, dc$idv1)
  fa <- survival::survfit(survival::Surv(time, status) ~ 1, data = d[d$treat == 1, ], weights = w)
  tt <- c(0, fa$time)
  expect_equal(step_at(v$x, v$y, tt), step_at(c(0, fa$time), c(1, fa$surv), tt), tolerance = 1e-10)
})

test_that("where S = 0 the band closes onto the curve", {
  d <- data.frame(time = c(1, 2, 3, 4, 1.5, 2.5, 3.5, 5), status = c(1, 1, 1, 1, 1, 0, 1, 0), treat = rep(0:1, each = 4))
  dc <- df_counting(d, "time", "status", "treat", get.cox = FALSE)
  v <- km_step_vertices(dc$at_points, dc$surv0, dc$idv0, se = sqrt(dc$sig2_surv0))
  k <- length(v$x)
  expect_equal(c(v$y[k], v$lower[k], v$upper[k]), c(0, 0, 0))
  expect_false(anyNA(c(v$lower, v$upper, v$band_y)))
})

test_that("plot_km_curves_counting() draws the helper's vertices", {
  d <- gbsg_er0()
  dc <- df_counting(d, "tte", "event", "treat", by.risk = 12)
  drawn <- list(); polys <- list()
  local_mocked_bindings(lines = function(x, y, ...) drawn[[length(drawn) + 1]] <<- list(x = x, y = y),
                        polygon = function(x, y, ...) polys[[length(polys) + 1]] <<- list(x = x, y = y))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  KM_plot_2sample_weighted_counting(dc, conf.int = TRUE, show.cox = FALSE, show.med = FALSE)
  v1 <- km_step_vertices(dc$at_points, dc$surv1, dc$idv1, se = sqrt(dc$sig2_surv1))
  v0 <- km_step_vertices(dc$at_points, dc$surv0, dc$idv0, se = sqrt(dc$sig2_surv0))
  expect_equal(drawn[[1]], list(x = v1$x, y = v1$y))
  expect_equal(drawn[[2]], list(x = v0$x, y = v0$y))
  expect_equal(polys[[1]], list(x = v0$band_x, y = v0$band_y))
  expect_equal(polys[[2]], list(x = v1$band_x, y = v1$band_y))
})

# B11: the two KM confidence bands are semi-transparent; their geometry is unchanged.

test_that("band fills have alpha below 1, keep their hues, and the vertices are the helper's", {
  g <- survival::gbsg[survival::gbsg$er <= 0, ]
  d <- data.frame(tte = g$rfstime / 30.4375, event = g$status, treat = g$hormon)
  dc <- df_counting(d, "tte", "event", "treat", by.risk = 12)
  polys <- list()
  local_mocked_bindings(polygon = function(x, y, col = NULL, ...) polys[[length(polys) + 1]] <<- list(x = x, y = y, col = col))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  KM_plot_2sample_weighted_counting(dc, conf.int = TRUE, show.cox = FALSE, show.med = FALSE)
  expect_length(polys, 2)
  rgba <- grDevices::col2rgb(vapply(polys, `[[`, "", "col"), alpha = TRUE)
  expect_true(all(rgba["alpha", ] < 255))
  expect_equal(rgba[1:3, 1], grDevices::col2rgb("lightgrey")[, 1])
  expect_equal(rgba[1:3, 2], grDevices::col2rgb("lightblue")[, 1])
  v0 <- km_step_vertices(dc$at_points, dc$surv0, dc$idv0, se = sqrt(dc$sig2_surv0))
  v1 <- km_step_vertices(dc$at_points, dc$surv1, dc$idv1, se = sqrt(dc$sig2_surv1))
  expect_equal(polys[[1]][c("x", "y")], list(x = v0$band_x, y = v0$band_y))
  expect_equal(polys[[2]][c("x", "y")], list(x = v1$band_x, y = v1$band_y))
})

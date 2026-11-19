test_that("V: df_counting() no longer defines a local check_km_curve()", {
  src <- paste(deparse(body(df_counting)), collapse = "\n")
  expect_false(grepl("check_km_curve <- function", src, fixed = TRUE))
  expect_identical(names(formals(check_km_curve)), c("S.KM", "group_name", "stop_on_error"))
})

test_that("V: df_counting() returns each element once", {
  r <- suppressWarnings(df_counting(ws_sim(100, 2), "time", "status", "treat", get.cox = FALSE))
  expect_false(anyDuplicated(names(r)) > 0)
  expect_true("strata" %in% names(r))
})

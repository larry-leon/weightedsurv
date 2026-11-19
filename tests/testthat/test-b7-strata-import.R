test_that("B7: strata() is imported from survival", {
  expect_true("strata" %in% getNamespaceImports("weightedsurv")$survival)
})

test_that("B7: stratified df_counting() runs", {
  v <- ws_veteran()
  r <- suppressWarnings(df_counting(v, "time", "status", "treat", strata.name = "celltype"))
  expect_true(is.finite(r$z.score_stratified))
})

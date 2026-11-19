# Shared, seeded test data (tests need only survival and testthat).
ws_sim <- function(n = 120, seed = 1, ties = FALSE, weighted = FALSE) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  d <- data.frame(treat = rep(0:1, length.out = n))
  d$time <- stats::rexp(n, 0.08 * exp(-0.3 * d$treat))
  cc <- stats::runif(n, 0, 30)
  d$status <- as.numeric(d$time <= cc)
  d$time <- pmin(d$time, cc)
  if (ties) d$time <- ceiling(d$time * 2) / 2
  d$w <- if (weighted) stats::runif(n, 0.3, 3) else 1
  d
}
ws_veteran <- function() {
  v <- survival::veteran
  v$treat <- v$trt - 1
  v
}

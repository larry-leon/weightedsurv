# W1: ps_weights()

test_that("W1: weights, IF and dw equal the reference, both types; the model route equals the formula route", {
  d <- ws_sim_ps(400, 11)
  for (type in c("ipw", "stabilized")) {
    pw <- ps_weights(treat ~ x1 + x2, data = d, type = type)
    ref <- ws_ps_ref(treat ~ x1 + x2, d, type)
    expect_s3_class(pw, "ps_weights")
    expect_identical(names(pw), c("weights", "type", "treat", "ps", "IF", "dw", "model", "call"))
    expect_lt(max(abs(pw$weights - ref$weights)), 1e-10)
    expect_lt(max(abs(pw$IF - ref$IF)), 1e-10)
    expect_lt(max(abs(pw$dw - ref$dw)), 1e-10)
    expect_identical(pw$treat, as.numeric(d$treat)); expect_identical(pw$type, type)
    expect_identical(dim(pw$IF), c(400L, 3L)); expect_identical(dim(pw$dw), c(400L, 3L))
    g <- stats::glm(treat ~ x1 + x2, family = stats::binomial(), data = d)
    pm <- ps_weights(model = g, type = type)
    expect_lt(max(abs(pm$weights - pw$weights)), 1e-10)
    expect_lt(max(abs(pm$IF - pw$IF)), 1e-10); expect_lt(max(abs(pm$dw - pw$dw)), 1e-10)
    expect_equal(pw$ps, as.numeric(stats::fitted(g)), tolerance = 1e-12)
  }
})

test_that("W1: dw and IF equal numerical derivatives", {
  d <- ws_sim_ps(400, 11); n <- nrow(d); A <- d$treat
  pw <- ps_weights(treat ~ x1 + x2, data = d)
  ctl <- stats::glm.control(epsilon = 1e-13, maxit = 100)
  g <- stats::glm(treat ~ x1 + x2, family = stats::binomial(), data = d, control = ctl); b <- stats::coef(g); xm <- stats::model.matrix(g)
  wof <- function(bb) { pp <- as.numeric(stats::plogis(xm %*% bb)); ifelse(A == 1, 1 / pp, 1 / (1 - pp)) }
  num <- sapply(seq_along(b), function(j) { h <- 1e-6; bp <- b; bp[j] <- b[j] + h; bm <- b; bm[j] <- b[j] - h; (wof(bp) - wof(bm)) / (2 * h) })
  expect_lt(max(abs(pw$dw - num)) / max(abs(num)), 1e-5)
  ii <- c(1, 7, 50, 200, 399)
  ifn <- t(sapply(ii, function(i) { wt <- rep(1, n); wt[i] <- 1 + 1e-4
    (stats::coef(stats::glm(treat ~ x1 + x2, family = stats::quasibinomial(), data = d, weights = wt, control = ctl)) - b) / 1e-4 }))
  expect_lt(max(abs(pw$IF[ii, ] - ifn)) / max(abs(ifn)), 1e-3)
})

test_that("W1: stabilised weights and dw are the ipw ones times the arm proportion; IF identical", {
  d <- ws_sim_ps(400, 11); p1 <- mean(d$treat); sc <- ifelse(d$treat == 1, p1, 1 - p1)
  pw <- ps_weights(treat ~ x1 + x2, data = d, type = "ipw"); pst <- ps_weights(treat ~ x1 + x2, data = d, type = "stabilized")
  expect_lt(max(abs(pst$weights - sc * pw$weights)), 1e-12)
  expect_lt(max(abs(pst$dw - sc * pw$dw)), 1e-12)
  expect_identical(pst$IF, pw$IF)
})

test_that("W1: aliased columns are dropped from IF and dw", {
  d <- ws_sim_ps(300, 12); d$x3 <- 2 * d$x1
  pw <- suppressWarnings(ps_weights(treat ~ x1 + x2 + x3, data = d))
  expect_identical(ncol(pw$IF), 3L); expect_identical(colnames(pw$dw), c("(Intercept)", "x1", "x2"))
  expect_lt(max(abs(pw$IF - ps_weights(treat ~ x1 + x2, data = d)$IF)), 1e-8)
})

test_that("W1: each misuse is an error that names the problem", {
  d <- ws_sim_ps(300, 13); n <- nrow(d)
  expect_error(ps_weights(model = stats::lm(treat ~ x1, data = d)), "must be a glm fit")
  expect_error(ps_weights(model = stats::glm(treat ~ x1, family = stats::binomial(link = "probit"), data = d)), "logistic regression.*probit")
  expect_error(ps_weights(model = stats::glm(time ~ x1, family = stats::gaussian(), data = d)), "logistic regression.*gaussian")
  expect_error(ps_weights(model = stats::glm(treat ~ x1, family = stats::binomial(), data = d, weights = rep(2, n))), "prior weights")
  expect_error(ps_weights(model = stats::glm(treat ~ x1 + offset(x2), family = stats::binomial(), data = d)), "offset")
  expect_error(ps_weights(model = suppressWarnings(stats::glm(treat ~ x1 + x2, family = stats::binomial(), data = d, control = stats::glm.control(maxit = 1)))), "did not converge")
  expect_error(ps_weights(treat ~ x1, data = transform(d, x1 = ifelse(seq_len(n) == 3, NA, x1))), "1 row\\(s\\) were dropped.*missing")
  dsep <- d; dsep$x9 <- ifelse(dsep$treat == 1, 1, -1) + 0.001 * dsep$x1
  expect_error(ps_weights(treat ~ x9, data = dsep), "did not converge|numerically 0 or 1")
  expect_error(ps_weights(treat ~ x1, data = transform(d, treat = 1)), "both values 0 and 1")
  expect_error(ps_weights(treat ~ x1, data = d, model = stats::glm(treat ~ x1, family = stats::binomial(), data = d)), "not both")
  expect_error(ps_weights(treat ~ x1), "both needed")
  expect_error(ps_weights(treat ~ x1, data = d, type = "overlap"), "should be one of")
  for (m in list(tryCatch(ps_weights(model = stats::lm(treat ~ x1, data = d)), error = conditionMessage),
                 tryCatch(ps_weights(treat ~ x1, data = transform(d, treat = 1)), error = conditionMessage)))
    expect_false(grepl("TRUE/FALSE|subscript out of bounds", m))
})

test_that("W1: print shows type, n, coefficients and by arm the size, weight range and effective sample size", {
  d <- ws_sim_ps(300, 13); pw <- ps_weights(treat ~ x1 + x2, data = d, type = "stabilized")
  out <- utils::capture.output(r <- print(pw))
  expect_identical(r, pw)
  expect_match(out[1], "stabilized.*n = 300.*3 coefficients")
  expect_true(any(grepl("ess", out)) && any(grepl("treat = 1", out)) && any(grepl("treat = 0", out)))
  w1 <- pw$weights[d$treat == 1]
  expect_true(any(grepl(format(sum(w1)^2 / sum(w1^2), digits = 4), out, fixed = TRUE)))
})

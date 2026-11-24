# The tolerances below are statistical, not numerical: they bound sampling error in an
# estimated derivative, not floating-point noise. They are tight enough that the
# bandwidth-scaled quantity the estimator used to return (roughly 1.4x the derivative
# here) fails them comfortably.

test_that("d = 1 estimates the derivative itself, not a bandwidth-scaled multiple", {
  d <- sim_known_erf_data(slope = 2)

  e1 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 41L, d = 1)

  # Interior anchors only; local polynomial derivative estimates are biased at the edges.
  interior <- e1$values[10:32]
  est <- e1$erf(interior, vcov = FALSE)$est

  expect_equal(mean(est), 2, tolerance = 0.05)
  expect_true(all(abs(est - 2) < 0.25))
})

test_that("d = 1 tracks a derivative that varies over the exposure", {
  d <- sim_known_erf_data(quad = TRUE)

  e1 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 41L, d = 1)

  interior <- e1$values[10:32]
  est <- e1$erf(interior, vcov = FALSE)$est

  # theta(a) = a^2, so theta'(a) = 2a.
  expect_equal(est, 2 * interior, tolerance = 0.1)
})

test_that("d = 1 agrees with the numerical derivative of the d = 0 curve", {
  d <- sim_known_erf_data(quad = TRUE)

  e0 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 41L, d = 0)
  e1 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 41L, d = 1)

  v <- e0$values
  est0 <- e0$erf(v, vcov = FALSE)$est

  # Central difference of the ERF, against the directly estimated derivative.
  num <- (est0[-(1:2)] - est0[-((length(v) - 1):length(v))]) / (v[-(1:2)] - v[-((length(v) - 1):length(v))])
  est1 <- e1$erf(v[-c(1, length(v))], vcov = FALSE)$est

  interior <- 8:30
  expect_equal(est1[interior], num[interior], tolerance = 0.15)
})

test_that("the local polynomial degree follows d", {
  d <- sim_known_erf_data(slope = 2)

  # p = d + 1, and the bandwidth is selected for that degree, so the two differ.
  e0 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, d = 0)
  e1 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, d = 1)

  expect_false(isTRUE(all.equal(e0$bw, e1$bw)))

  # get_G() builds the design matrix, so its width is the degree actually used.
  expect_identical(ncol(get_G(rnorm(10), 0L + 1L)), 2L)
  expect_identical(ncol(get_G(rnorm(10), 1L + 1L)), 3L)
})

test_that("d is validated and stored", {
  expect_identical(process_d(0), 0L)
  expect_identical(process_d(1), 1L)

  for (bad in list(2, -1, 0.5, 10)) {
    expect_error(process_d(bad), regexp = "must be", ignore.case = TRUE)
  }

  d <- sim_known_erf_data(n = 400)
  expect_identical(estimate_erf(Y ~ A | X1 + X2, data = d, n = 11L, d = 1)$d, 1L)
  expect_identical(estimate_erf(Y ~ A | X1 + X2, data = d, n = 11L)$d, 0L)

  expect_error(estimate_erf(Y ~ A | X1 + X2, data = d, n = 11L, d = 2),
               regexp = "must be", ignore.case = TRUE)
})

test_that("standard errors are on the same scale as the derivative", {
  d <- sim_known_erf_data(slope = 2)

  e1 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, d = 1)

  a <- e1$values[11]
  r <- e1$erf(a, vcov = TRUE)
  se <- sqrt(drop(r$vcov))

  # A Wald interval for the derivative should cover the truth. If the estimate were
  # rescaled without its standard error (or the reverse), this would fail.
  expect_lt(abs(r$est - 2) / se, 4)
})

test_that("a tree grown from a derivative ERF uses the derivative at its nodes", {
  d <- sim_known_erf_data(n = 800, slope = 2)
  d$M <- factor(rbinom(nrow(d), 1, 0.5))

  e1 <- estimate_erf(Y ~ A | X1 + X2 + M, data = d, n = 21L, d = 1)
  tr <- rpart_erf(e1, ~ M + X1, max_depth = 2, min_samples = 100)

  # The root node covers the whole sample, so its ERF is the one passed in.
  expect_equal(tr[[1]]$erf(e1$values, vcov = FALSE)$est,
               e1$erf(e1$values, vcov = FALSE)$est)

  # And it is the derivative, near 2, not the ERF, which rises across the range.
  expect_equal(mean(tr[[1]]$erf(e1$values[8:14], vcov = FALSE)$est), 2,
               tolerance = 0.2)
})

test_that("a derivative tree splits on the slope modifier, not the level modifier", {
  # Growing a tree from a `d = 1` ERF applies the chosen objective to the derivative
  # curves rather than to the curves, which is what makes a derivative-based objective
  # available without adding one. `Lv` shifts the level of the curve and `Sl` changes its
  # slope, so the two orders should select different modifiers.
  set.seed(21)

  n <- 900
  X1 <- rnorm(n)
  Lv <- factor(rbinom(n, 1, 0.5))
  Sl <- factor(rbinom(n, 1, 0.5))
  A <- 1 + 0.5 * X1 + rnorm(n)
  Y <- 10 * (Lv == "1") + A * (1 + 2 * (Sl == "1")) + X1 + rnorm(n)
  d <- data.frame(Y, A, X1, Lv, Sl)

  root_var <- function(e, metric) {
    tr <- rpart_erf(e, ~ Lv + Sl, metric = metric, max_depth = 1, min_samples = 150)
    tr$tree$split_var
  }

  e0 <- estimate_erf(Y ~ A | X1 + Lv + Sl, data = d, n = 21L, d = 0)
  e1 <- estimate_erf(Y ~ A | X1 + Lv + Sl, data = d, n = 21L, d = 1)

  expect_identical(root_var(e0, "L2"), "Lv")
  expect_identical(root_var(e0, "AUC"), "Lv")

  expect_identical(root_var(e1, "L2"), "Sl")
  expect_identical(root_var(e1, "AUC"), "Sl")

  # And the node curves under `d = 1` are the subgroup slopes, 1 and 3.
  tr <- rpart_erf(e1, ~ Lv + Sl, metric = "L2", max_depth = 1, min_samples = 150)
  interior <- e1$values[7:15]
  slopes <- sort(vapply(2:3, function(i) mean(tr[[i]]$erf(interior, vcov = FALSE)$est),
                        numeric(1L)))

  expect_equal(slopes, c(1, 3), tolerance = 0.3)
})

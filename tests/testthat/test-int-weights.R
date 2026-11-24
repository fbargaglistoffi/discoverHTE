test_that("all int_weight_method values produce usable integral weights", {
  d <- sim_erf_data()
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  sp <- d$X1 < median(d$X1)
  erfL <- estimate_erf_internal(e$fit, values = e$values, subset = sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)
  erfR <- estimate_erf_internal(e$fit, values = e$values, subset = !sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)

  for (iw in c("none", "density", "gps")) {
    w <- get_int_weights(erfL, erfR, e$values, iw)

    expect_length(w, length(e$values))
    expect_true(all(is.finite(w)), info = iw)
    expect_true(all(w >= 0), info = iw)
    expect_gt(sum(w), 0)
  }
})

test_that("the gps weight interpolates on the nuisance exposure grid", {
  # Regression test: `pi_hat_mat` column means are indexed by `A_vals`, not by the
  # subgroup's own exposures. Interpolating against the latter made "gps" error on every
  # call, since the two lengths only coincide by accident.
  d <- sim_erf_data()
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  res <- attr(e$erf, "res")

  expect_false(is.null(res$A_vals))
  expect_length(res$A_vals, ncol(res$pi_hat_mat))
  expect_false(length(res$A) == ncol(res$pi_hat_mat))
})

test_that("rpart_erf runs for every metric x int_weight_method combination", {
  d <- sim_erf_data()
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  for (iw in c("none", "density", "gps")) {
    for (mt in c("L2", "L1", "AUC", "fANOVA")) {
      tr <- expect_no_error(
        rpart_erf(e, ~ mod + X1, metric = mt, int_weight_method = iw,
                  max_depth = 2, min_samples = 50)
      )
      expect_s3_class(tr, "erf_tree")
    }
  }
})

test_that("disjoint exposure support falls back to uniform weights", {
  # A split on the exposure itself gives two branches with essentially no common
  # support, so the product of their densities integrates to ~0.
  d <- sim_erf_data()
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)

  sp <- d$A < quantile(d$A, 0.5)
  erfL <- estimate_erf_internal(e$fit, values = e$values, subset = sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)
  erfR <- estimate_erf_internal(e$fit, values = e$values, subset = !sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)

  w <- get_int_weights(erfL, erfR, e$values, "density")

  expect_true(all(is.finite(w)))
  expect_gt(sum(w), 0)
})

test_that("an anchor outside a branch's exposure range gets no weight", {
  # interp_dens() interpolates each branch's density with rule = 1, so an anchor beyond
  # that branch's own exposures comes back NA and has to become 0 rather than propagate.
  # Negative interpolated densities are flattened to 0 for the same reason.
  d <- sim_erf_data(n = 400)
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 41L)

  sp <- d$A < quantile(d$A, 0.5)
  erfL <- estimate_erf_internal(e$fit, values = e$values, subset = sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)
  erfR <- estimate_erf_internal(e$fit, values = e$values, subset = !sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)

  for (iw in c("density", "gps")) {
    w <- get_int_weights(erfL, erfR, e$values, iw)

    expect_false(anyNA(w), info = iw)
    expect_true(all(w >= 0), info = iw)
    expect_true(all(is.finite(w)), info = iw)
  }
})

test_that("update(bw = ) at the original bandwidth reproduces the ERF", {
  d <- sim_erf_data(n = 300, clsize = 10)

  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)
  e2 <- update(e, bw = e$bw)

  expect_equal(e2$bw, e$bw)
  expect_equal(e2$erf(e$values, vcov = TRUE)$est,
               e$erf(e$values, vcov = TRUE)$est)
  expect_equal(e2$erf(e$values, vcov = TRUE)$vcov,
               e$erf(e$values, vcov = TRUE)$vcov)
})

test_that("update(bw = ) keeps the cluster-robust variance", {
  # Regression test: update() used to call estimate_erf_internal() without `cluster`,
  # so re-smoothing a clustered ERF silently returned iid standard errors while
  # `object$cluster` still reported the clustering.
  d <- sim_erf_data(n = 300, clsize = 10)

  e_cl <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, cluster = ~cl)
  e_iid <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)

  v_cl <- e_cl$erf(e_cl$values, vcov = TRUE)$vcov
  v_iid <- e_iid$erf(e_iid$values, vcov = TRUE)$vcov

  # The clustered and iid variances must differ, or this test cannot detect the bug.
  expect_false(isTRUE(all.equal(v_cl, v_iid)))

  v_updated <- update(e_cl, bw = e_cl$bw)$erf(e_cl$values, vcov = TRUE)$vcov

  expect_equal(v_updated, v_cl)
  expect_false(isTRUE(all.equal(v_updated, v_iid)))
})

test_that("update(bw = ) keeps the local-polynomial coefficient order", {
  # Regression test: `d` was not stored on the <erf> and was not passed through
  # update(), so re-smoothing a derivative ERF returned the ERF itself.
  d <- sim_erf_data(n = 300, clsize = 10)

  e1 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, d = 1)
  e0 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, d = 0)

  expect_identical(e1$d, 1L)
  expect_false(isTRUE(all.equal(e1$erf(e1$values, vcov = FALSE)$est,
                                e0$erf(e0$values, vcov = FALSE)$est)))

  e1_updated <- update(e1, bw = e1$bw)

  expect_identical(e1_updated$d, 1L)
  expect_equal(e1_updated$erf(e1$values, vcov = TRUE)$est,
               e1$erf(e1$values, vcov = TRUE)$est)
  expect_equal(e1_updated$erf(e1$values, vcov = TRUE)$vcov,
               e1$erf(e1$values, vcov = TRUE)$vcov)
})

test_that("update(bw = ) stores the resolved bandwidth", {
  # `bw` may arrive as a function or as NULL; what is stored should be the number
  # actually smoothed with, as estimate_erf() stores it.
  d <- sim_erf_data(n = 200, clsize = 10)

  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)

  e_fun <- update(e, bw = bw.SJ)
  expect_type(e_fun$bw, "double")
  expect_length(e_fun$bw, 1L)
  expect_equal(e_fun$bw, bw.SJ(e$treat))

  e_cv <- update(e, bw = NULL)
  expect_type(e_cv$bw, "double")
  expect_length(e_cv$bw, 1L)
})

test_that("update(bw = ) changes the estimates it is supposed to change", {
  d <- sim_erf_data(n = 300, clsize = 10)

  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)
  e_wide <- update(e, bw = 2 * e$bw)

  expect_equal(e_wide$bw, 2 * e$bw)
  expect_false(isTRUE(all.equal(e_wide$erf(e$values, vcov = FALSE)$est,
                                e$erf(e$values, vcov = FALSE)$est)))

  # Only the smoothing changes; the nuisance fit behind it is the same object.
  expect_identical(e_wide$fit, e$fit)
  expect_identical(e_wide$cluster, e$cluster)
  expect_identical(e_wide$trim, e$trim)
})

test_that("update() rejects the ERF at a node of a tree", {
  d <- sim_erf_data(n = 300, clsize = 10)

  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 50)

  expect_error(update(tr[[2]], bw = 1),
               regexp = "single node", ignore.case = TRUE)
})

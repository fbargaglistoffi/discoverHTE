test_that("fANOVA weights each child's deviation by its share of the node", {
  d <- sim_erf_data(n = 400)
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)

  # Deliberately lopsided: 10% vs 90%.
  sp <- d$X1 < quantile(d$X1, 0.1)

  erfL <- estimate_erf_internal(e$fit, values = e$values, subset = sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)
  erfR <- estimate_erf_internal(e$fit, values = e$values, subset = !sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)

  iw <- get_int_weights(erfL, erfR, e$values, "density")

  got <- get_metric(erfL, erfR, values = e$values, metric = "fANOVA",
                    int_weights = iw, erf_parent = e$erf)

  # Reconstruct both the weighted and the (former) unweighted objective.
  integrate_erf <- function(x, y, w) {
    wa <- (w[-1L] + w[-length(w)]) / 2
    sum(diff(x) * (y[-1L] + y[-length(y)]) * wa) / 2
  }

  parent <- e$erf(e$values, vcov = FALSE)$est
  L2_L <- integrate_erf(e$values, (erfL(e$values, vcov = FALSE)$est - parent)^2, iw)
  L2_R <- integrate_erf(e$values, (erfR(e$values, vcov = FALSE)$est - parent)^2, iw)

  nL <- sum(attr(erfL, "res")[["weights"]])
  nR <- sum(attr(erfR, "res")[["weights"]])

  expect_equal(got, (nL * L2_L + nR * L2_R) / (nL + nR))

  # The weighting is what keeps a small, wildly-deviating child from scoring like a
  # large one, so it must be strictly below the old unweighted sum here.
  expect_lt(got, L2_L + L2_R)
})

test_that("fANOVA weights sum to one across the two children", {
  # The split is exhaustive, so weighting by child share must reduce to a plain average
  # when both children deviate from the parent by the same amount.
  d <- sim_erf_data(n = 400)
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)

  sp <- d$X1 < median(d$X1)

  erfL <- estimate_erf_internal(e$fit, values = e$values, subset = sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)
  erfR <- estimate_erf_internal(e$fit, values = e$values, subset = !sp,
                                bw = e$bw, vcov = FALSE, verbose = FALSE)

  nL <- sum(attr(erfL, "res")[["weights"]])
  nR <- sum(attr(erfR, "res")[["weights"]])

  expect_equal(nL + nR, nrow(d))
})

test_that("a candidate split that cannot be scored does not abort the tree", {
  d <- sim_erf_data()
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  # get_metric() is called once per candidate; make it fail for some of them. The tree
  # must still be built from whatever candidates remain scoreable.
  calls <- 0L
  local_mocked_bindings(
    get_metric = function(...) {
      calls <<- calls + 1L
      if (calls %% 2L == 0L) stop("synthetic failure")
      1 / calls
    }
  )

  tr <- expect_no_error(rpart_erf(e, ~ mod + X1, max_depth = 1, min_samples = 50))
  expect_s3_class(tr, "erf_tree")
  expect_gt(calls, 1L)
})

test_that("a node whose every candidate fails becomes a leaf", {
  d <- sim_erf_data()
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  local_mocked_bindings(
    get_metric = function(...) stop("synthetic failure")
  )

  tr <- expect_no_error(rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 50))

  expect_identical(tr$tree$type, "leaf")
  expect_identical(tr$tree$reason, "no_valid_modifers")
})

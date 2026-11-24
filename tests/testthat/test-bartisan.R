# Keep the sampler small; these tests check the plumbing, not the fit quality.
bart_method <- function(...) {
  method_bartisan(num_burn = 50L, num_draws = 100L, num_trees = 20L, ...)
}

test_that("method_bartisan recovers a linear ERF", {
  skip_if_not_installed("bartisan")

  d <- sim_erf_data(n = 400)

  set.seed(42)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, method = bart_method(),
                    n_grid = 30L, n = 25L)

  expect_s3_class(e, "erf")
  expect_identical(attr(e$method, "method_name"), "bartisan")

  grid <- seq(-1.5, 3, length.out = 10)
  est <- e$erf(grid, vcov = FALSE)$est

  expect_true(all(is.finite(est)))

  # Y = A * (1 + 1.5 * (mod == "1")) + ..., with mod balanced, so the ERF is
  # linear with slope 1.75.
  expect_equal(unname(coef(lm(est ~ grid))[2L]), 1.75, tolerance = 0.15)
})

test_that("the nuisance components from method_bartisan are well formed", {
  skip_if_not_installed("bartisan")

  d <- sim_erf_data(n = 200)

  set.seed(42)
  nuis <- bart_method()(Y ~ A | X1 + X2 + mod, data = d, n_grid = 15L)

  expect_named(nuis, c("Y", "A", "A_vals", "pidens", "mupred", "weights", "cluster"))
  expect_length(nuis$A_vals, 15L)

  pd <- nuis$pidens(seq_len(200L))
  mp <- nuis$mupred(seq_len(200L))

  expect_length(pd$den_dens, 200L)
  expect_identical(dim(pd$pi_hat_mat), c(200L, 15L))
  expect_length(mp$mu_hat, 200L)
  expect_identical(dim(mp$mu_hat_mat), c(200L, 15L))

  # A density is nonnegative, and the one at the observed exposure is the
  # diagonal of the grid matrix in spirit, so it must be positive somewhere.
  expect_true(all(pd$den_dens >= 0))
  expect_true(all(pd$pi_hat_mat >= 0))
  expect_true(any(pd$den_dens > 0))
})

test_that("method_bartisan works with cross-fitting and with weights", {
  skip_if_not_installed("bartisan")

  d <- sim_erf_data(n = 300)
  grid <- seq(-1, 2.5, length.out = 6)

  set.seed(7)
  e_cf <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, method = bart_method(),
                       n_grid = 20L, n = 15L, crossfit = 3L)

  expect_true(all(is.finite(e_cf$erf(grid, vcov = FALSE)$est)))

  w <- runif(nrow(d), 0.5, 2)

  set.seed(7)
  e_unw <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, method = bart_method(),
                        n_grid = 20L, n = 15L)
  set.seed(7)
  e_w <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, weights = w,
                      method = bart_method(), n_grid = 20L, n = 15L)

  expect_equal(e_w$weights, w, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(e_unw$erf(grid, vcov = FALSE)$est,
                                e_w$erf(grid, vcov = FALSE)$est)))
})

test_that("families supplied to method_bartisan are used", {
  skip_if_not_installed("bartisan")

  d <- sim_erf_data(n = 200)
  grid <- seq(-1, 2.5, length.out = 6)

  set.seed(42)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d,
                    method = bart_method(outcome = list(family = gaussian()),
                                         exposure = list(family = gaussian())),
                    n_grid = 15L, n = 15L)

  expect_true(all(is.finite(e$erf(grid, vcov = FALSE)$est)))

  # A name that is neither a `bartisan()` nor a `bartisan_control()` argument is
  # an error rather than being silently dropped.
  expect_error(estimate_erf(Y ~ A | X1 + X2 + mod, data = d,
                            method = bart_method(not_an_argument = 1),
                            n_grid = 15L, n = 15L))
})

test_that("an erf fit with method_bartisan can be split by rpart_erf", {
  skip_if_not_installed("bartisan")

  d <- sim_erf_data(n = 400)

  set.seed(42)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, method = bart_method(),
                    n_grid = 20L, n = 21L)

  tr <- rpart_erf(e, ~ mod + X1 + X2, max_depth = 1, min_samples = 60)

  expect_s3_class(tr, "erf_tree")
  expect_identical(tr$tree$type, "internal")
  expect_identical(tr$tree$split_var, "mod")
  expect_identical(attr(tr$method, "method_name"), "bartisan")
})

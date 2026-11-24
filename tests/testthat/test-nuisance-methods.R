# Every nuisance method targets the same estimand, so each should recover a known ERF on
# data where the outcome and exposure models are both correctly specified. The tolerances
# are statistical: they bound sampling error at these sample sizes, not floating-point
# noise. An unadjusted regression of Y on A has a slope near 2.3 on this data, so a method
# that failed to adjust for the confounders would not pass.

test_that("method_glm recovers a known linear ERF", {
  d <- sim_known_erf_data(n = 1000, slope = 2)

  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L, method = method_glm())

  interior <- e$values[6:16]
  expect_equal(e$erf(interior, vcov = FALSE)$est, 2 * interior, tolerance = 0.1)
  expect_identical(attr(e$method, "method_name"), "glm")
})

test_that("method_glm recovers a nonlinear ERF", {
  d <- sim_known_erf_data(n = 1500, quad = TRUE)

  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)

  # method_glm() puts a natural cubic spline on the exposure, so a quadratic ERF is
  # within its outcome model.
  interior <- e$values[6:16]
  expect_equal(e$erf(interior, vcov = FALSE)$est, interior^2, tolerance = 0.2)
})

test_that("method_glm passes families to the model it names", {
  d <- sim_known_erf_data(n = 600)
  set.seed(3)
  d$count <- rpois(nrow(d), exp(0.3 * d$A))

  # `exposure` defaults to `outcome`, so a Poisson outcome family would otherwise be
  # handed to a continuous exposure.
  e <- estimate_erf(count ~ A | X1 + X2, data = d, n = 15L,
                    method = method_glm(outcome = list(family = poisson()),
                                        exposure = list(family = gaussian())))

  expect_s3_class(e, "erf")
  expect_true(all(is.finite(e$erf(e$values, vcov = FALSE)$est)))

  # The ERF of a count outcome rising in exposure must itself rise.
  est <- e$erf(e$values[4:12], vcov = FALSE)$est
  expect_gt(est[length(est)], est[1L])
})

test_that("method_gamlss recovers a known linear ERF", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  d <- sim_known_erf_data(n = 600, slope = 2)

  e <- suppressWarnings(
    estimate_erf(Y ~ A | X1 + X2, data = d, n = 15L,
                 method = method_gamlss(outcome = list(family = gamlss.dist::NO()),
                                        exposure = list(family = gamlss.dist::NO())))
  )

  interior <- e$values[4:12]
  expect_equal(e$erf(interior, vcov = FALSE)$est, 2 * interior, tolerance = 0.2)
  expect_identical(attr(e$method, "method_name"), "gamlss")
})

# A binary outcome takes the binomial family with a logit link, so the outcome model is a
# logistic regression with a penalized spline on the exposure. The data are generated
# from a logistic model linear in the exposure, so the spline should be penalized to a
# line and the fit should match glm().
sim_binary_erf_data <- function(n = 600, seed = 21) {
  set.seed(seed)

  X1 <- rnorm(n)
  X2 <- rnorm(n)
  A <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)
  Y <- rbinom(n, 1, plogis(-1 + 0.6 * A + 0.5 * X1 - 0.4 * X2))

  data.frame(Y, A, X1, X2)
}

test_that("method_gamlss fits a logistic regression to a binary outcome", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  d <- sim_binary_erf_data()

  # The default outcome family for a two-valued outcome, so no family is supplied for it.
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 11L,
                    method = method_gamlss(exposure = list(family = gamlss.dist::NO())))

  mu <- e$fit$mupred(seq_len(nrow(d)))

  # The fitted probabilities are those of a logistic regression on the same terms.
  g <- glm(Y ~ A + X1 + X2, data = d, family = binomial())
  expect_equal(mu$mu_hat, unname(fitted(g)), tolerance = 0.01)

  # Every prediction, including those across the exposure grid, is a probability.
  expect_true(all(mu$mu_hat > 0 & mu$mu_hat < 1))
  expect_true(all(mu$mu_hat_mat > 0 & mu$mu_hat_mat < 1))

  # The ERF is on the probability scale, rises with the exposure as the outcome model
  # does, and carries finite standard errors.
  est <- e$erf(e$values, vcov = TRUE)
  expect_true(all(est$est > 0 & est$est < 1))
  expect_gt(est$est[9L], est$est[3L])
  expect_true(all(is.finite(diag(est$vcov)) & diag(est$vcov) > 0))
})

test_that("method_gamlss recovers a known ERF for a binary outcome", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  d <- sim_binary_erf_data(n = 1500)

  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 11L, range = .9,
                    method = method_gamlss(exposure = list(family = gamlss.dist::NO())))

  # The true ERF is the mean of the true conditional probability over the covariates.
  truth <- vapply(e$values, function(a) {
    mean(plogis(-1 + 0.6 * a + 0.5 * d$X1 - 0.4 * d$X2))
  }, numeric(1L))

  expect_equal(e$erf(e$values, vcov = FALSE)$est, truth, tolerance = 0.1)
})

test_that("method_gamlss cross-fits a binary outcome", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  d <- sim_binary_erf_data()

  set.seed(4)
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 9L, crossfit = 2,
                    method = method_gamlss(exposure = list(family = gamlss.dist::NO())))

  mu <- e$fit$mupred(seq_len(nrow(d)))
  expect_true(all(mu$mu_hat > 0 & mu$mu_hat < 1))
  expect_true(all(is.finite(e$erf(e$values, vcov = FALSE)$est)))
})

test_that("method_gamlss accepts a family as an object, a function, or a name", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  d <- sim_binary_erf_data(n = 400)

  fit_with <- function(outcome_family, exposure_family) {
    estimate_erf(Y ~ A | X1 + X2, data = d, n = 9L,
                 method = method_gamlss(outcome = list(family = outcome_family),
                                        exposure = list(family = exposure_family)))
  }

  e_obj <- fit_with(gamlss.dist::BI(), gamlss.dist::NO())
  e_fun <- fit_with(gamlss.dist::BI, gamlss.dist::NO)
  e_chr <- fit_with("BI", "NO")

  expect_equal(e_fun$erf(e_fun$values, vcov = FALSE)$est,
               e_obj$erf(e_obj$values, vcov = FALSE)$est)
  expect_equal(e_chr$erf(e_chr$values, vcov = FALSE)$est,
               e_obj$erf(e_obj$values, vcov = FALSE)$est)
})

test_that("method_gamlss rejects a family from stats", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  d <- sim_binary_erf_data(n = 200)

  expect_error(estimate_erf(Y ~ A | X1 + X2, data = d, n = 9L,
                            method = method_gamlss(outcome = list(family = binomial()))),
               "BI()", fixed = TRUE)
})

test_that("method_gamlss works with a one-parameter outcome family", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  # The Poisson has a mean parameter only, so its fitted parameters come back as a single
  # column, the case that a binomial outcome on the exposure grid also produces.
  d <- sim_known_erf_data(n = 500)
  set.seed(5)
  d$count <- rpois(nrow(d), exp(0.2 + 0.3 * d$A))

  e <- estimate_erf(count ~ A | X1 + X2, data = d, n = 9L,
                    method = method_gamlss(outcome = list(family = gamlss.dist::PO()),
                                           exposure = list(family = gamlss.dist::NO())))

  est <- e$erf(e$values, vcov = FALSE)$est
  expect_true(all(is.finite(est) & est > 0))
  expect_gt(est[7L], est[3L])
})

# On this data, an ST3 outcome model fits `tau` at or below 1 for some units at high
# exposures, where ST3 has no mean. It is the case that made the ST3 default fail.
sim_heavy_tau_data <- function() {
  set.seed(1)
  n <- 2000

  X1 <- rnorm(n)
  X2 <- rnorm(n)
  A <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)
  Yb <- rbinom(n, 1, plogis(-1 + 0.6 * A + 0.5 * X1 - 0.4 * X2))
  Y <- 2 * A + X1 + X2 + rnorm(n)

  data.frame(Y, A, X1, X2)[1:600, ]
}

test_that("the default continuous outcome family always has a mean", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  d <- sim_heavy_tau_data()

  # SST, whose degrees of freedom are kept above 2, so every unit has a finite mean.
  e <- suppressWarnings(
    estimate_erf(Y ~ A | X1 + X2, data = d, n = 9L, method = method_gamlss())
  )

  mu <- e$fit$mupred(seq_len(nrow(d)))
  expect_true(all(is.finite(mu$mu_hat)))
  expect_true(all(is.finite(mu$mu_hat_mat)))

  interior <- e$values[3:7]
  expect_equal(e$erf(interior, vcov = FALSE)$est, 2 * interior, tolerance = 0.1)
})

test_that("an outcome family with no finite mean is reported by name", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  # Whether the fitted `tau` crosses 1 depends on the optimizer's path, which can differ
  # across BLAS libraries, so this is not run on CRAN.
  skip_on_cran()

  d <- sim_heavy_tau_data()

  expect_error(
    suppressWarnings(
      estimate_erf(Y ~ A | X1 + X2, data = d, n = 9L,
                   method = method_gamlss(outcome = list(family = gamlss.dist::ST3())))
    ),
    "no finite mean"
  )
})

test_that("method_superlearner recovers a known linear ERF", {
  skip_if_not_installed("SuperLearner")

  d <- sim_known_erf_data(n = 600, slope = 2)

  # A seed is required: the cross-validation splits are random.
  set.seed(11)
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 15L,
                    method = method_superlearner(
                      SL.library = c("SL.mean", "SL.glm")))

  interior <- e$values[4:12]
  expect_equal(e$erf(interior, vcov = FALSE)$est, 2 * interior, tolerance = 0.25)
  expect_identical(attr(e$method, "method_name"), "superlearner")
})

test_that("the nuisance components are well formed for every available method", {
  d <- sim_known_erf_data(n = 400, slope = 2)

  methods <- list(glm = method_glm())

  if (rlang::is_installed("SuperLearner")) {
    set.seed(11)
    methods$superlearner <- method_superlearner(SL.library = c("SL.mean", "SL.glm"))
  }

  for (nm in names(methods)) {
    e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 11L, method = methods[[nm]])
    res <- attr(e$erf, "res")

    # The outcome-mean and density matrices are both n x n_grid, and the density is a
    # density, so it cannot be negative.
    expect_identical(dim(res$mu_hat_mat), dim(res$pi_hat_mat), info = nm)
    expect_identical(nrow(res$mu_hat_mat), length(res$A), info = nm)
    expect_true(all(is.finite(res$mu_hat_mat)), info = nm)
    expect_true(all(res$pi_hat_mat >= 0), info = nm)
    expect_true(all(is.finite(res$pseudo_Y)), info = nm)
  }
})

test_that("trim bounds the density ratio without changing the ERF much", {
  d <- sim_known_erf_data(n = 800, slope = 2)

  e_trim <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 15L, trim = 0.99)
  e_none <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 15L, trim = 1)

  expect_identical(e_trim$trim, 0.99)
  expect_identical(e_none$trim, 1)

  # Trimming trades a little bias for variance, so the two must differ but agree broadly.
  est_trim <- e_trim$erf(e_trim$values[4:12], vcov = FALSE)$est
  est_none <- e_none$erf(e_none$values[4:12], vcov = FALSE)$est

  expect_false(isTRUE(all.equal(est_trim, est_none)))
  expect_equal(est_trim, est_none, tolerance = 0.25)

  expect_error(estimate_erf(Y ~ A | X1 + X2, data = d, n = 11L, trim = 1.5))
  expect_error(estimate_erf(Y ~ A | X1 + X2, data = d, n = 11L, trim = 0))
})

test_that("sampling weights are carried into the fit", {
  d <- sim_known_erf_data(n = 600, slope = 2)
  set.seed(5)
  d$w <- runif(nrow(d), 0.5, 2)

  e_w <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 15L, weights = d$w)
  e_u <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 15L)

  expect_equal(e_w$weights, d$w)
  expect_true(all(e_u$weights == 1))

  # A string naming a column in `data` is accepted too, and gives the same fit.
  e_s <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 15L, weights = "w")
  expect_equal(e_s$weights, e_w$weights)
  expect_equal(e_s$erf(e_s$values, vcov = FALSE)$est,
               e_w$erf(e_w$values, vcov = FALSE)$est)

  expect_error(estimate_erf(Y ~ A | X1 + X2, data = d, n = 11L, weights = ~w))

  # Weighting changes the estimate, but both still recover the ERF.
  interior <- e_w$values[4:12]
  expect_false(isTRUE(all.equal(e_w$erf(interior, vcov = FALSE)$est,
                                e_u$erf(interior, vcov = FALSE)$est)))
  expect_equal(e_w$erf(interior, vcov = FALSE)$est, 2 * interior, tolerance = 0.25)
})

test_that("method_superlearner cross-fits every part of the exposure density", {
  # Regression test: the kernel density of standardized residuals is the last piece of
  # the fitted exposure density, and it used to be estimated from the *evaluation* fold's
  # residuals while the two SuperLearner models were trained on the training fold. That
  # put the evaluation fold inside the nuisance estimate used to evaluate it.
  #
  # The leak is not observable through the public API, so what is pinned here is the
  # boundary condition that makes the fix safe: with no cross-fitting the training and
  # evaluation rows are the same, so the density is unchanged and results must not move.
  skip_if_not_installed("SuperLearner")

  d <- sim_known_erf_data(n = 400, slope = 2)
  lib <- c("SL.mean", "SL.glm")

  fit <- function(cf) {
    set.seed(11)
    estimate_erf(Y ~ A | X1 + X2, data = d, n = 15L, crossfit = cf,
                 method = method_superlearner(SL.library = lib))
  }

  e_none <- fit(0)
  e_false <- fit(FALSE)
  e_one <- fit(1)

  expect_identical(e_false$erf(e_none$values, vcov = FALSE)$est,
                   e_none$erf(e_none$values, vcov = FALSE)$est)
  expect_identical(e_one$erf(e_none$values, vcov = FALSE)$est,
                   e_none$erf(e_none$values, vcov = FALSE)$est)
})

test_that("the cross-fitted SuperLearner exposure density stays a density", {
  skip_if_not_installed("SuperLearner")

  d <- sim_known_erf_data(n = 400, slope = 2)
  d$id <- factor(rep(seq_len(40), each = 10))

  set.seed(11)
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 15L, crossfit = 3, cluster = ~id,
                    method = method_superlearner(SL.library = c("SL.mean", "SL.glm")))

  res <- attr(e$erf, "res")

  # Every unit gets a finite, positive conditional density at its own exposure, and the
  # grid of conditional densities carries no gaps. Fitting the residual density on a
  # fold that does not cover the evaluation points would show up here as NAs.
  expect_false(anyNA(res$pi_hat_mat))
  expect_true(all(is.finite(res$pi_hat_mat)))
  expect_true(all(res$pi_hat_mat >= 0))

  # Each row of pi_hat_mat is pi(a | X_i) over the exposure grid, so it integrates to
  # about 1 over the range the grid spans.
  grid <- res$A_vals
  integrate_row <- function(y) {
    sum(diff(grid) * (y[-1L] + y[-length(y)])) / 2
  }
  mass <- apply(res$pi_hat_mat, 1L, integrate_row)

  expect_equal(mean(mass), 1, tolerance = 0.1)
  expect_true(all(mass > 0.5 & mass < 1.5))

  # The derivative of the ERF is still recovered, so the corrected density has not
  # broken the estimator it feeds.
  expect_equal(mean(e$erf(e$values[4:12], vcov = FALSE)$est), 2 * mean(e$values[4:12]),
               tolerance = 0.3)
})

#' Nuisance Functional Estimation
#' @name erf_methods
#'
#' @description
#' The output of these functions is supplied to the `method` argument of [estimate_erf()] to determine how the nuisance functionals are estimated when estimating the ERF.
#'
#' @param outcome a list of arguments passed to the outcome model. See Details for more information.
#' @param exposure a list of arguments passed to the exposure model. See Details for more information.
#' @param ... arguments passed to both the exposure and outcome models.
#'
#' @return
#' An `<erf_nuisance>` object, which is a function that produces the nuisance functional estimates required for estimating the pseudo-outcomes and influence functions for the ERF.
#'
#' @details
#' The following sections describe what each nuisance estimator does. The required components are estimates of the conditional outcome mean given the exposure and covariates and estimates of the conditional density of exposure given the covariates. From these components, the other components required for estimation of the pseudo-outcomes and influence functions can be computed directly.
#'
#' ## `method_glm`
#'
#' Estimates the conditional outcome mean and the conditional density of exposure with generalized linear models. First, a model for the exposure given the covariates is fit to arrive at estimates for the parameters of the conditional exposure distribution. These are used directly to compute the conditional exposure density based on the model family. Second, a model for the outcome given the exposure and covariates is fit to arrive at estimates for the conditional outcome mean.
#'
#' For the exposure, a linear model of the covariates is fit. For the outcome, a linear model with a natural cubic spline term for the exposure with 4 df (using \pkgfun{splines}{ns}) is fit.
#'
#' ## `method_superlearner`
#'
#' Estimates the conditional outcome mean and the conditional density of exposure with SuperLearner using functions in the \CRANpkg{SuperLearner} package, which stacks multiple machine learning model predictions together to optimize cross-validated prediction error. First, a model for the exposure given the covariates is fit to arrive at estimates of the conditional exposure mean. Second, a model for the log of squared exposure residuals given the covariates is fit to arrive estimates of the conditional exposure variance. Both of these models are fit by calling \pkgfun{SuperLearner}{SuperLearner} with the arguments supplied to `exposure`. The conditional exposure density is estimated by applying kernel density estimation to the ratio of the conditional exposure estimates to the square root of the conditional variance estimates. Finally, the conditional outcome mean is estimated by calling \pkgfun{SuperLearner}{SuperLearner} with the arguments supplied to `outcome`.
#'
#' The `SL.library` arguments are the most important to specify. This controls which learners are used to train the exposure and outcome models. Ideally, this should contain a mix of parametric and nonparametric models. See \pkgfun{SuperLearner}{listWrappers} for options. The default is
#' \preformatted{SL.library = c("SL.gam", "SL.glm", "SL.glm.interaction",
#'                "SL.mean", "SL.rpart")}
#' for both the exposure and outcome models.
#'
#' Note that due to randomness in how the cross-validation splits are performed, a seed must be set with [set.seed()] to ensure replicability.
#'
#' ## `method_gamlss`
#'
#' Estimates the conditional outcome mean and the conditional density of exposure with generalized additive models for location, scale, and shape (GAMLSS) using functions in the \CRANpkg{gamlss} package. GAMLSS are related to generalized linear models but include models for the whole conditional density of the response variable, not just the mean. They can also include additive smoothing terms for continuous covariates. First, a model for the exposure given the covariates is fit to arrive at estimates for the parameters of the conditional exposure distribution. These are used directly to compute the conditional exposure density. Second, a model for the outcome given the exposure and covariates is fit to arrive at estimates for the parameters of the conditional outcome distribution, from which the conditional outcome mean is estimated.
#'
#' For each model, a distribution must be chosen. Ideally, this should be one that matches the nature of the response variables, e.g., a model with natural number support when the response is a count. See \pkgfun{gamlss.dist}{gamlss.family} for options; a family can be supplied as a family object (e.g., `gamlss.dist::NO()`), the function that creates one (`gamlss.dist::NO`), or its name (`"NO"`). Families from \pkg{stats}, such as `binomial()`, are not accepted. For the outcome, the default is \pkgfun{gamlss.dist}{BI}, the binomial distribution with a logit link, when the outcome takes two values, and \pkgfun{gamlss.dist}{SST}, the skew-t type 3 distribution with 4 parameters and support for the real numbers, otherwise. A binary outcome must be coded as 0 and 1, and its outcome model is then a logistic regression. `SST` is parameterized by its mean and standard deviation and restricts its degrees of freedom to be greater than 2, so the conditional outcome mean always exists; \pkgfun{gamlss.dist}{ST3}, the same distribution in its original parameterization, has no mean where its fitted degrees of freedom are 1 or less, and an error is thrown when the outcome model yields such a fit. For the exposure, the default is \pkgfun{gamlss.dist}{ST3}, for which only the density is needed.
#'
#' For the exposure, a linear model of the covariates is fit for each parameter in the response distribution. For the outcome, a linear model with a penalized spline term for the exposure (using \pkgfun{gamlss}{pb}) is fit for each parameter in the response distribution.
#'
#' ## `method_bartisan`
#'
#' Estimates the conditional outcome mean and the conditional density of exposure with Bayesian additive regression trees (BART) using functions in the \CRANpkg{bartisan} package. BART is a sum-of-trees model in which a regularizing prior keeps each tree a weak learner, so nonlinearity and interactions are found rather than specified. First, a model for the exposure given the covariates is fit, and the conditional exposure density is computed directly from the fitted response distribution at each exposure value. Second, a model for the outcome given the exposure and covariates is fit to arrive at estimates for the conditional outcome mean.
#'
#' For each model, a response distribution must be chosen; see `help("bartisan-families", package = "bartisan")` for options. For the outcome, the default is `binomial()` when the outcome takes two values and `gaussian()` otherwise. For the exposure, the default is \pkgfun{bartisan}{gaussian_ls}, the location-scale Gaussian, in which a second forest models the log of the conditional standard deviation so that the estimated density can vary in spread as well as in location across the covariates. When the exposure is bounded or discrete, a family matching its support should be supplied instead. Note that `exposure` does not default to `outcome` for this method, so a family supplied to `outcome` applies to the outcome model alone and the exposure keeps its own default.
#'
#' Arguments supplied to `...` are passed to both models, which is the convenient way to reach the sampler and prior settings documented at \pkgfun{bartisan}{bartisan_control}, as in `method_bartisan(num_trees = 20, num_draws = 200)`. Note that both models are predicted at every combination of unit and exposure grid value, so the cost of this method grows with `n_grid`; lowering `n_grid` in the call to [estimate_erf()] is the first thing to try when a fit takes too long.
#'
#' Note that due to randomness in the posterior sampling, a seed must be set with [set.seed()] to ensure replicability.
#'
#' ## `method_negbinomcens`
#'
#' Estimates the nuisance functionals for a count outcome measured over multiple observations clustered within `unit`s (e.g., age bands within a county-year), returning them on an age-standardized rate scale at the unit level. This method is intended for age-standardized mortality/incidence rates where the raw counts may be interval-censored (e.g., suppressed small cells). It requires \CRANpkg{mgcv}.
#'
#' The conditional outcome mean is estimated with a negative-binomial generalized additive model (via \pkgfun{mgcv}{bam}) for the count outcome, with a thin-plate/cubic-regression spline term for the exposure, spline terms for continuous covariates, parametric terms for discrete covariates and `age`, and a `log(population)` offset. Interval-censored cells (a missing outcome with finite `censor_lo`/`censor_hi` bounds) are handled with an EM loop that imputes each censored cell's conditional expectation `E[y | censor_lo <= y <= censor_hi]` under the fitted negative-binomial. Predicted and observed counts are converted to rates (per 100,000), age-standardized within each `unit` using `std_weights`, and it is these unit-level quantities that enter the pseudo-outcome. The g-computation curve is obtained by predicting the count model across the exposure grid, which makes no assumption about the structure of the outcome model. (A faster shortcut exploiting the log-additive structure of this particular model is retained as commented-out code in the source.)
#'
#' The conditional density of exposure is estimated as in `method_glm()` (a generalized linear model of the exposure on the covariates, with the density implied by the model family), but fit on the *unit-level* dataset (one row per `unit`) since the exposure is constant within a unit.
#'
#' Clustering is supplied to [estimate_erf()] rather than to this method, as `estimate_erf(..., cluster = ~county)`; it must not vary within a `unit`. It is aggregated to the unit level and used both to form cross-fitting folds over whole clusters and to compute cluster-robust standard errors.
#'
#' @seealso [estimate_erf()]
#'
#' @examples
#' set.seed(1234)
#' n <- 500
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' A <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)
#' Y <- 2 * A + X1 + X2 + rnorm(n)
#' d <- data.frame(Y, A, X1, X2)
#'
#' # Generalized linear models for both nuisance functionals, the default
#' e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 25L,
#'                   method = method_glm())
#' e
#'
#' # A count outcome. `exposure` defaults to `outcome`, so the exposure model
#' # is given a family of its own rather than inheriting the Poisson one.
#' d$count <- rpois(n, exp(0.3 * A))
#'
#' e_pois <- estimate_erf(count ~ A | X1 + X2, data = d, n = 25L,
#'                        method = method_glm(outcome = list(family = poisson()),
#'                                            exposure = list(family = gaussian())))
#' e_pois
#'
#' @examplesIf rlang::is_installed("SuperLearner")
#' # Stacked learners for both models; a seed is needed for the
#' # cross-validation splits
#' set.seed(1234)
#' e_sl <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 25L,
#'                      method = method_superlearner(
#'                        SL.library = c("SL.mean", "SL.glm")))
#' e_sl
#'
#' @export
#' @rdname erf_methods
method_glm <- function(outcome = list(),
                       exposure = list(),
                       ...) {
  arg::arg_list(outcome)
  arg::arg_list(exposure)

  out <- function(formula, data, weights = NULL, crossfit = 0L, n_grid = 100L,
                  cluster = NULL) {
    .c(Y, A, X) %=% process_formula(formula, data) #multi-assign

    A_vals <- seq(min(A), max(A), length.out = n_grid)

    n <- length(A)

    if (is_null(weights)) {
      weights <- rep.int(1.0, n)
    }

    if (...length() > 0) {
      dots <- list(...)

      exposure <- utils::modifyList(exposure, dots)
      outcome <- utils::modifyList(outcome, dots)
    }

    if (is_null(exposure[["family"]])) {
      exposure$family <- gaussian()
    }

    if (is_null(outcome[["family"]])) {
      outcome$family <- if (length(unique(Y)) <= 2L) binomial() else gaussian()
    }

    dat <- cbind(Y = Y, A = A, X)

    # estimate nuisance functions via glm
    exposure_formula <- formula(dat[-1L])

    outcome_formula <- formula(dat[-2L]) |>
      update(as.formula(sprintf(". ~ . + splines::ns(%s, df = 4)", names(dat)[2L])))

    fit_models <- function(id_train, id_eval = id_train) {
      ## Train models on id_train
      exposure <- exposure |>
        utils::modifyList(list(formula = exposure_formula,
                               data = dat,
                               subset = id_train,
                               weights = weights))

      pimod <- do.call(stats::glm, exposure,
                       quote = TRUE)

      make_dens_fun <- function(pimod) {
        fam <- family(pimod)$family
        if (fam == "Gamma") {
          alpha_hat <- 1 / summary(pimod)$dispersion
          function(a, fitted) dgamma(a, shape = alpha_hat, rate = alpha_hat / fitted)
        }
        else if (fam == "poisson") {
          function(a, fitted) dpois(a, fitted)
        }
        else {
          sd_hat <- sigma(pimod)
          function(a, fitted) dnorm(a, fitted, sd_hat)
        }
      }

      dens_fun <- make_dens_fun(pimod)

      outcome <- outcome |>
        utils::modifyList(list(formula = outcome_formula,
                               data = dat,
                               subset = id_train,
                               weights = weights))

      mumod <- do.call(stats::glm, outcome,
                       quote = TRUE)

      ## Compute predictions on id_eval
      n_eval <- length(id_eval)

      # Pred units
      A_grid_eval <- rep(A_vals, each = n_eval)
      X_grid_eval <- ss(X, id_eval) |> ss(rep(seq_len(n_eval), n_grid))

      AX_grid_eval <- cbind(A = A_grid_eval, X_grid_eval) |>
        as.data.frame()

      den_dens <- predict(pimod, newdata = ss(dat, id_eval), type = "response") |>
        dens_fun(a = ss(A, id_eval))

      pi_hat_mat <- predict(pimod, newdata = AX_grid_eval, type = "response") |>
        dens_fun(a = A_grid_eval) |>
        matrix(nrow = n_eval)

      mu_hat <- predict(mumod, newdata = ss(dat, id_eval), type = "response")

      mu_hat_mat <- predict(mumod, newdata = AX_grid_eval, type = "response") |>
        matrix(nrow = n_eval)

      list(
        den_dens = den_dens,
        pi_hat_mat = pi_hat_mat,
        mu_hat = mu_hat,
        mu_hat_mat = mu_hat_mat
      )
    }

    if (crossfit <= 1L) {
      fit_out <- fit_models(seq_len(n))

      den_dens <- fit_out$den_dens
      pi_hat_mat <- fit_out$pi_hat_mat
      mu_hat <- fit_out$mu_hat
      mu_hat_mat <- fit_out$mu_hat_mat
    }
    else {
      den_dens <- numeric(n)
      pi_hat_mat <- matrix(0, nrow = n, ncol = n_grid)
      mu_hat <- numeric(n)
      mu_hat_mat <- matrix(0, nrow = n, ncol = n_grid)

      folds <- make_folds(crossfit, n, cluster[[1L]])

      for (i in seq_len(crossfit)) {
        id_train <- folds %!=% i
        id_eval <-  folds %==% i

        fit_out <- fit_models(id_train, id_eval)

        den_dens[id_eval] <- fit_out$den_dens
        pi_hat_mat[id_eval, ] <- fit_out$pi_hat_mat
        mu_hat[id_eval] <- fit_out$mu_hat
        mu_hat_mat[id_eval, ] <- fit_out$mu_hat_mat
      }
    }

    pidens <- function(subset) {
      list(
        den_dens = ss(den_dens, subset),
        pi_hat_mat = ss(pi_hat_mat, subset)
      )
    }

    mupred <- function(subset) {
      list(
        ## Predicted Y for each unit
        mu_hat = ss(mu_hat, subset),

        ## Matrix of predicted Y for each X * A_grid grid unit
        mu_hat_mat = ss(mu_hat_mat, subset)
      )
    }

    list(Y = Y,
         A = A,
         A_vals = A_vals,
         pidens = pidens,
         mupred = mupred,
         weights = weights,
         # Already aligned to the analysis rows, since this method does not aggregate.
         cluster = cluster)
  }

  attr(out, "method_name") <- "glm"
  class(out) <- c(class(out), "erf_nuisance")

  out
}

#' @export
#' @rdname erf_methods
method_superlearner <- function(outcome = list(SL.library = c("SL.gam", "SL.glm", "SL.glm.interaction",
                                                              "SL.mean", "SL.rpart")),
                                exposure = list(SL.library = c("SL.gam", "SL.glm", "SL.glm.interaction",
                                                               "SL.mean", "SL.rpart")),
                                ...) {
  rlang::check_installed("SuperLearner",
                         "to use `method_superlearner()`.")

  arg::arg_list(outcome)
  arg::arg_list(exposure)

  out <- function(formula, data, weights = NULL, crossfit = 0L, n_grid = 100L,
                  cluster = NULL) {
    .c(Y, A, X) %=% process_formula(formula, data) #multi-assign

    A_vals <- seq(min(A), max(A), length.out = n_grid)

    n <- length(A)

    if (is_null(weights)) {
      weights <- rep.int(1.0, n)
    }

    if (...length() > 0) {
      dots <- list(...)

      exposure <- utils::modifyList(exposure, dots)
      outcome <- utils::modifyList(outcome, dots)
    }

    out_family <- if (length(unique(Y)) <= 2L) binomial() else gaussian()

    AX <- cbind(A = A, X) |>
      as.data.frame()

    fit_models <- function(id_train, id_eval = id_train) {
      n_eval <- length(id_eval)

      # Pred units (need to define here to include in SuperLearner call)
      A_grid_eval <- rep(A_vals, each = n_eval)
      X_grid_eval <- ss(X, id_eval) |> ss(rep(seq_len(n_eval), n_grid))

      AX_grid_eval <- cbind(A = A_grid_eval, X_grid_eval) |>
        as.data.frame()

      ## Train models on id_train
      args <- exposure |>
        utils::modifyList(list(Y = A[id_train], X = ss(X, id_train),
                               newX = rbind(X, X_grid_eval),
                               obsWeights = ss(weights, id_train),
                               family = gaussian(),
                               env = asNamespace("SuperLearner")))

      pimod <- do.call(SuperLearner::SuperLearner, args)


      args$Y <- log((A[id_train] - pimod$SL.predict[1:n][id_train])^2)
      pi2mod <- do.call(SuperLearner::SuperLearner, args)

      args <- outcome |>
        utils::modifyList(list(Y = Y[id_train], X = ss(AX, id_train),
                               newX = rbind(ss(AX, id_eval), AX_grid_eval),
                               obsWeights = ss(weights, id_train),
                               family = out_family,
                               env = asNamespace("SuperLearner")))

      mumod <- do.call(SuperLearner::SuperLearner, args)

      ## Compute predictions on id_eval
      pimod.vals_i <- pimod$SL.predict[1:n][id_eval]
      pimod.vals_grid <- pimod$SL.predict[-(1:n)]

      pi2mod.vals_i <- pi2mod$SL.predict[1:n][id_eval] |> exp() |> sqrt()
      pi2mod.vals_grid <- pi2mod$SL.predict[-(1:n)] |> exp() |> sqrt()

      A_std_eval <- (A[id_eval] - pimod.vals_i) / pi2mod.vals_i
      A_std_grid_eval <- (A_grid_eval - pimod.vals_grid) / pi2mod.vals_grid

      # The density of the standardized residuals is the last piece of the fitted
      # exposure density, so it is estimated on the training units, like the two models
      # above it. Fitting it on `id_eval` would put the evaluation fold inside the
      # nuisance estimate used to evaluate that same fold, which is what cross-fitting
      # exists to prevent. Without cross-fitting `id_train` and `id_eval` are the same
      # rows and this is the same density as before.
      pimod.vals_train <- pimod$SL.predict[1:n][id_train]
      pi2mod.vals_train <- pi2mod$SL.predict[1:n][id_train] |> exp() |> sqrt()

      A_std_train <- (A[id_train] - pimod.vals_train) / pi2mod.vals_train

      # `A_std_eval` always lies inside the range of `A_std_grid_eval`: the grid spans
      # `A_vals` at each evaluation unit's own fitted mean and scale, and that unit's own
      # exposure is one of the values in between. So the interpolator built over the grid
      # range covers every point the density is evaluated at.
      A_std_dens <- density(A_std_train, bw = "SJ",
                            weights = weights[id_train] / sum(weights[id_train]),
                            from = min(A_std_grid_eval),
                            to = max(A_std_grid_eval)) |>
        with(approxfun(x = x, y = y))

      den_dens <- A_std_dens(A_std_eval) / pi2mod.vals_i

      pi_hat_mat <- matrix(A_std_dens(A_std_grid_eval) / pi2mod.vals_grid,
                           nrow = n_eval)

      mu_hat <- mumod$SL.predict[1:n_eval]

      mu_hat_mat <- matrix(mumod$SL.predict[-(1:n_eval)],
                           nrow = n_eval)

      list(
        den_dens = den_dens,
        pi_hat_mat = pi_hat_mat,
        mu_hat = mu_hat,
        mu_hat_mat = mu_hat_mat
      )
    }

    if (crossfit <= 1L) {
      fit_out <- fit_models(seq_len(n))

      den_dens <- fit_out$den_dens
      pi_hat_mat <- fit_out$pi_hat_mat
      mu_hat <- fit_out$mu_hat
      mu_hat_mat <- fit_out$mu_hat_mat
    }
    else {
      den_dens <- numeric(n)
      pi_hat_mat <- matrix(0, nrow = n, ncol = n_grid)
      mu_hat <- numeric(n)
      mu_hat_mat <- matrix(0, nrow = n, ncol = n_grid)

      folds <- make_folds(crossfit, n, cluster[[1L]])

      for (i in seq_len(crossfit)) {
        id_train <- folds %!=% i
        id_eval <-  folds %==% i

        fit_out <- fit_models(id_train, id_eval)

        den_dens[id_eval] <- fit_out$den_dens
        pi_hat_mat[id_eval, ] <- fit_out$pi_hat_mat
        mu_hat[id_eval] <- fit_out$mu_hat
        mu_hat_mat[id_eval, ] <- fit_out$mu_hat_mat
      }
    }

    pidens <- function(subset) {
      list(
        den_dens = ss(den_dens, subset),
        pi_hat_mat = ss(pi_hat_mat, subset)
      )
    }

    mupred <- function(subset) {
      list(
        ## Predicted Y for each unit
        mu_hat = ss(mu_hat, subset),

        ## Matrix of predicted Y for each X * A_grid grid unit
        mu_hat_mat = ss(mu_hat_mat, subset)
      )
    }

    list(Y = Y,
         A = A,
         A_vals = A_vals,
         pidens = pidens,
         mupred = mupred,
         weights = weights,
         # Already aligned to the analysis rows, since this method does not aggregate.
         cluster = cluster)
  }

  attr(out, "method_name") <- "superlearner"
  class(out) <- c(class(out), "erf_nuisance")

  out
}

#' @export
#' @rdname erf_methods
method_gamlss <- function(outcome = list(),
                          exposure = list(),
                          ...) {
  rlang::check_installed(c("gamlss", "gamlss.dist"),
                         "to use `method_gamlss()`.")

  arg::arg_list(outcome)
  arg::arg_list(exposure)

  out <- function(formula, data, weights = NULL, crossfit = 0L, n_grid = 100L,
                  cluster = NULL) {
    .c(Y, A, X) %=% process_formula(formula, data) #multi-assign

    A_vals <- seq(min(A), max(A), length.out = n_grid)

    n <- length(A)

    if (is_null(weights)) {
      weights <- rep.int(1.0, n)
    }

    if (...length() > 0) {
      dots <- list(...)

      exposure <- utils::modifyList(exposure, dots)
      outcome <- utils::modifyList(outcome, dots)
    }

    # The families are stored as gamlss family objects, not the functions that make them,
    # because the formulas below read which parameters each family has from
    # `family$parameters`. A binary outcome gets the binomial with a logit link, so its
    # mean model is a logistic regression; stats::binomial() is not a gamlss family.
    #
    # Any other outcome gets SST, the skew t type 3 distribution parameterized by its mean
    # and standard deviation. Its `tau` link keeps the degrees of freedom above 2, so the
    # mean always exists. ST3, the same distribution in its original parameterization,
    # has no mean where the fitted `tau` is at or below 1, and with `tau` modeled on the
    # covariates and the exposure that happens for some units. The exposure model needs
    # only a density, which ST3 has for every `tau`.
    if (is_null(exposure[["family"]])) {
      exposure$family <- gamlss.dist::ST3()
    }

    if (is_null(outcome[["family"]])) {
      outcome$family <- if (length(unique(Y)) <= 2L) gamlss.dist::BI() else gamlss.dist::SST()
    }

    exposure$family <- as_gamlss_family(exposure$family, "exposure")
    outcome$family <- as_gamlss_family(outcome$family, "outcome")

    # estimate nuisance functions via gamlss
    dat <- cbind(Y = Y, A = A, X)

    mu.formula_exposure <- formula(dat[-1L])

    sigma.formula_exposure <- {
      if (isTRUE(exposure$family$parameters$sigma)) reformulate(names(dat)[-(1:2)])
      else ~1
    }

    nu.formula_exposure <- {
      if (isTRUE(exposure$family$parameters$nu)) reformulate(names(dat)[-(1:2)])
      else ~1
    }

    tau.formula_exposure <- {
      if (isTRUE(exposure$family$parameters$tau)) reformulate(names(dat)[-(1:2)])
      else ~1
    }


    mu.formula_outcome <- formula(dat) |>
      update(as.formula(sprintf(". ~ . + gamlss::pb(%s)", names(dat)[2L])))

    sigma.formula_outcome <- {
      if (isTRUE(outcome$family$parameters$sigma)) reformulate(names(dat)[-1L]) |>
        update(as.formula(sprintf("~ . + gamlss::pb(%s)", names(dat)[2L])))
      else ~1
    }

    nu.formula_outcome <- {
      if (isTRUE(outcome$family$parameters$nu)) reformulate(names(dat)[-1L]) |>
        update(as.formula(sprintf("~ . + gamlss::pb(%s)", names(dat)[2L])))
      else ~1
    }

    tau.formula_outcome <- {
      if (isTRUE(outcome$family$parameters$tau)) reformulate(names(dat)[-1L]) |>
        update(as.formula(sprintf("~ . + gamlss::pb(%s)", names(dat)[2L])))
      else ~1
    }

    fit_models <- function(id_train, id_eval = id_train) {
      ## Train models on id_train
      exposure <- exposure |>
        utils::modifyList(list(formula = mu.formula_exposure,
                               sigma.formula = sigma.formula_exposure,
                               nu.formula = nu.formula_exposure,
                               tau.formula = tau.formula_exposure,
                               data = ss(dat, id_train),
                               weights = weights[id_train],
                               trace = FALSE))

      pimod <- do.call(gamlss::gamlss, exposure,
                       quote = TRUE)

      ### Compute conditional density from estimated params
      dens_mod <- get(paste0("d", pimod$family[1L]), envir = asNamespace("gamlss.dist"))
      dens_args <- pimod$parameters

      # `drop = FALSE` keeps a one-parameter family's coefficients a list (a data frame
      # row) rather than a bare number, which `do.call()` would reject.
      dens <- function(a, coefs) {
        vapply(seq_along(a), function(i) {
          do.call(dens_mod, c(list(x = a[i]), as.list(coefs[i, dens_args, drop = FALSE])))
        }, numeric(1L))
      }

      outcome <- outcome |>
        utils::modifyList(list(formula = mu.formula_outcome,
                               sigma.formula = sigma.formula_outcome,
                               nu.formula = nu.formula_outcome,
                               tau.formula = tau.formula_outcome,
                               data = ss(dat, id_train),
                               weights = weights[id_train],
                               trace = FALSE))

      mumod <- do.call(gamlss::gamlss, outcome,
                       quote = TRUE)

      ### Compute conditional mean from estimated params
      mean_fun <- gamlss.dist::gamlss.family(family(mumod)[1L])[["mean"]]
      mean_args <- rlang::fn_fmls_names(mean_fun)

      get_mu <- function(coefs) {
        # The binomial families take the number of trials, `bd`, which predictAll()
        # returns only when `newdata` contains the outcome. Where it does not (the
        # exposure grid), it is 1, as it is for a binary outcome.
        if ("bd" %in% mean_args && !"bd" %in% names(coefs)) {
          coefs[["bd"]] <- 1
        }

        coefs <- ss(coefs, j = intersect(names(coefs), mean_args))

        # As in `dens()`, `drop = FALSE` keeps a one-column row a list.
        vapply(seq_row(coefs), function(i) {
          do.call(mean_fun, as.list(coefs[i, , drop = FALSE]))
        }, numeric(1L))
      }

      ## Compute predictions on id_eval
      n_eval <- length(id_eval)

      # Pred units
      A_grid_eval <- rep(A_vals, each = n_eval)
      X_grid_eval <- ss(X, id_eval) |> ss(rep(seq_len(n_eval), n_grid))

      AX_grid_eval <- cbind(A = A_grid_eval, X_grid_eval) |>
        as.data.frame()

      den_dens <- gamlss::predictAll(pimod, output = "data.frame",
                                     newdata = ss(dat, id_eval),
                                     data = ss(dat, id_train)) |>
        dens(a = A[id_eval])

      pi_hat_mat <- gamlss::predictAll(pimod, output = "data.frame",
                                       newdata = AX_grid_eval,
                                       data = ss(dat, id_train)) |>
        dens(a = A_grid_eval) |>
        matrix(nrow = n_eval)

      mu_hat <- gamlss::predictAll(mumod, output = "data.frame",
                                   newdata = ss(dat, id_eval),
                                   data = ss(dat, id_train)) |>
        get_mu()

      mu_hat_mat <- gamlss::predictAll(mumod, output = "data.frame",
                                       newdata = AX_grid_eval,
                                       data = ss(dat, id_train)) |>
        get_mu() |>
        matrix(nrow = n_eval)

      # A family's mean can be undefined at some fitted parameter values (e.g., ST3 with
      # `tau` at or below 1), which would otherwise surface as a non-finite pseudo-outcome
      # deep in the local linear fit.
      if (!all(is.finite(mu_hat)) || !all(is.finite(mu_hat_mat))) {
        fam <- mumod$family[1L]

        arg::err(c("The outcome model's fitted {.val {fam}} distribution has no finite mean for some units.",
                   i = "Choose a family whose mean is always defined, such as {.code gamlss.dist::SST()} in place of {.code gamlss.dist::ST3()}."))
      }

      list(
        den_dens = den_dens,
        pi_hat_mat = pi_hat_mat,
        mu_hat = mu_hat,
        mu_hat_mat = mu_hat_mat
      )
    }

    if (crossfit <= 1L) {
      fit_out <- fit_models(seq_len(n))

      den_dens <- fit_out$den_dens
      pi_hat_mat <- fit_out$pi_hat_mat
      mu_hat <- fit_out$mu_hat
      mu_hat_mat <- fit_out$mu_hat_mat
    }
    else {
      den_dens <- numeric(n)
      pi_hat_mat <- matrix(0, nrow = n, ncol = n_grid)
      mu_hat <- numeric(n)
      mu_hat_mat <- matrix(0, nrow = n, ncol = n_grid)

      folds <- make_folds(crossfit, n, cluster[[1L]])

      for (i in seq_len(crossfit)) {
        id_train <- folds %!=% i
        id_eval <-  folds %==% i

        fit_out <- fit_models(id_train, id_eval)

        den_dens[id_eval] <- fit_out$den_dens
        pi_hat_mat[id_eval, ] <- fit_out$pi_hat_mat
        mu_hat[id_eval] <- fit_out$mu_hat
        mu_hat_mat[id_eval, ] <- fit_out$mu_hat_mat
      }
    }

    pidens <- function(subset) {
      list(
        den_dens = ss(den_dens, subset),
        pi_hat_mat = ss(pi_hat_mat, subset)
      )
    }

    mupred <- function(subset) {
      list(
        ## Predicted Y for each unit
        mu_hat = ss(mu_hat, subset),

        ## Matrix of predicted Y for each X * A_grid grid unit
        mu_hat_mat = ss(mu_hat_mat, subset)
      )
    }

    list(Y = Y,
         A = A,
         A_vals = A_vals,
         pidens = pidens,
         mupred = mupred,
         weights = weights,
         # Already aligned to the analysis rows, since this method does not aggregate.
         cluster = cluster)
  }

  attr(out, "method_name") <- "gamlss"
  class(out) <- c(class(out), "erf_nuisance")

  out
}

# A family for gamlss() as a gamlss family object. The function that makes one (e.g.,
# `gamlss.dist::BI`) or its name (`"BI"`) is converted, since method_gamlss() reads the
# family's parameters before fitting. A family from `stats`, such as `binomial()`, is not
# a gamlss family, and gamlss() rejects it with an uninformative error, so it is caught
# here.
as_gamlss_family <- function(family, model) {
  if (inherits(family, "family") && !inherits(family, "gamlss.family")) {
    arg::err(c("The {model} family must be a {.pkg gamlss.dist} family, not {.fn {family$family}} from {.pkg stats}.",
               i = "For a binary {model}, use {.code gamlss.dist::BI()}, the binomial with a logit link."))
  }

  if (is.function(family) || rlang::is_string(family)) {
    family <- gamlss.dist::as.gamlss.family(family)
  }

  if (!inherits(family, "gamlss.family")) {
    arg::err("the {model} family must be a {.pkg gamlss.dist} family, such as {.code gamlss.dist::NO()}")
  }

  family
}

#' @export
#' @rdname erf_methods
method_bartisan <- function(outcome = list(),
                            exposure = list(),
                            ...) {
  rlang::check_installed("bartisan", "to use `method_bartisan()`.")

  arg::arg_list(outcome)
  arg::arg_list(exposure)

  out <- function(formula, data, weights = NULL, crossfit = 0L, n_grid = 100L,
                  cluster = NULL) {
    .c(Y, A, X) %=% process_formula(formula, data) #multi-assign

    A_vals <- seq(min(A), max(A), length.out = n_grid)

    n <- length(A)

    if (is_null(weights)) {
      weights <- rep.int(1.0, n)
    }

    if (...length() > 0) {
      dots <- list(...)

      exposure <- utils::modifyList(exposure, dots)
      outcome <- utils::modifyList(outcome, dots)
    }

    # Location-scale Gaussian for the exposure so the density can vary in spread
    # as well as in location; binary outcomes get a binomial, the rest Gaussian.
    if (is_null(exposure[["family"]])) {
      exposure$family <- bartisan::gaussian_ls()
    }

    if (is_null(outcome[["family"]])) {
      outcome$family <- if (length(unique(Y)) <= 2L) binomial() else gaussian()
    }

    dat <- cbind(Y = Y, A = A, X)

    # estimate nuisance functions via BART; no spline term for the exposure
    # because the trees supply the nonlinearity themselves
    exposure_formula <- formula(dat[-1L])

    outcome_formula <- formula(dat)

    fit_models <- function(id_train, id_eval = id_train) {
      ## Train models on id_train
      exposure <- exposure |>
        utils::modifyList(list(formula = exposure_formula,
                               data = ss(dat, id_train),
                               weights = ss(weights, id_train)))

      pimod <- do.call(bartisan::bartisan, exposure,
                       quote = TRUE)

      outcome <- outcome |>
        utils::modifyList(list(formula = outcome_formula,
                               data = ss(dat, id_train),
                               weights = ss(weights, id_train)))

      mumod <- do.call(bartisan::bartisan, outcome,
                       quote = TRUE)

      ## Compute predictions on id_eval
      n_eval <- length(id_eval)

      # Pred units
      A_grid_eval <- rep(A_vals, each = n_eval)
      X_grid_eval <- ss(X, id_eval) |> ss(rep(seq_len(n_eval), n_grid))

      AX_grid_eval <- cbind(A = A_grid_eval, X_grid_eval) |>
        as.data.frame()

      # The conditional exposure density is read off the fitted model at the
      # exposure value carried by `newdata`, so no residual-based KDE is needed.
      den_dens <- predict(pimod, newdata = ss(dat, id_eval),
                          type = "density")

      pi_hat_mat <- predict(pimod, newdata = AX_grid_eval,
                            type = "density") |>
        matrix(nrow = n_eval)

      mu_hat <- predict(mumod, newdata = ss(dat, id_eval),
                        type = "response")

      mu_hat_mat <- predict(mumod, newdata = AX_grid_eval,
                            type = "response") |>
        matrix(nrow = n_eval)

      list(
        den_dens = den_dens,
        pi_hat_mat = pi_hat_mat,
        mu_hat = mu_hat,
        mu_hat_mat = mu_hat_mat
      )
    }

    if (crossfit <= 1L) {
      fit_out <- fit_models(seq_len(n))

      den_dens <- fit_out$den_dens
      pi_hat_mat <- fit_out$pi_hat_mat
      mu_hat <- fit_out$mu_hat
      mu_hat_mat <- fit_out$mu_hat_mat
    }
    else {
      den_dens <- numeric(n)
      pi_hat_mat <- matrix(0, nrow = n, ncol = n_grid)
      mu_hat <- numeric(n)
      mu_hat_mat <- matrix(0, nrow = n, ncol = n_grid)

      folds <- make_folds(crossfit, n, cluster[[1L]])

      for (i in seq_len(crossfit)) {
        id_train <- folds %!=% i
        id_eval <-  folds %==% i

        fit_out <- fit_models(id_train, id_eval)

        den_dens[id_eval] <- fit_out$den_dens
        pi_hat_mat[id_eval, ] <- fit_out$pi_hat_mat
        mu_hat[id_eval] <- fit_out$mu_hat
        mu_hat_mat[id_eval, ] <- fit_out$mu_hat_mat
      }
    }

    pidens <- function(subset) {
      list(
        den_dens = ss(den_dens, subset),
        pi_hat_mat = ss(pi_hat_mat, subset)
      )
    }

    mupred <- function(subset) {
      list(
        ## Predicted Y for each unit
        mu_hat = ss(mu_hat, subset),

        ## Matrix of predicted Y for each X * A_grid grid unit
        mu_hat_mat = ss(mu_hat_mat, subset)
      )
    }

    list(Y = Y,
         A = A,
         A_vals = A_vals,
         pidens = pidens,
         mupred = mupred,
         weights = weights,
         # Already aligned to the analysis rows, since this method does not aggregate.
         cluster = cluster)
  }

  attr(out, "method_name") <- "bartisan"
  class(out) <- c(class(out), "erf_nuisance")

  out
}

#' @param age,population for `method_negbinomcens()`, the names of the variables in
#'   `data` giving the age band (modeled as a factor in the count outcome model) and
#'   the population denominator (entered as a `log()` offset).
#' @param unit for `method_negbinomcens()`, the name of the variable in `data`
#'   identifying the unit at which the ERF is estimated (e.g., a county-year).
#'   Observations are age-standardized and aggregated to units after the outcome
#'   model is fit; the ERF and its (cluster-robust) variance are computed at the unit
#'   level.
#' @param censor_lo,censor_hi for `method_negbinomcens()`, names of variables in
#'   `data` giving the lower and upper bounds for interval-censored outcome cells. A
#'   row with a missing outcome but finite bounds is treated as censored and its
#'   contribution imputed via EM.
#' @param std_weights for `method_negbinomcens()`, an optional named vector of
#'   age-standardization weights; defaults to the 2000 US standard population for age
#'   bands `"65-74"`, `"75-84"`, and `"85+"`, and to equal weights otherwise.
#' @param k_exposure,k_conf for `method_negbinomcens()`, the basis dimensions for the
#'   exposure and continuous-covariate smooths in the count outcome model.
#' @param em_iter for `method_negbinomcens()`, the maximum number of EM iterations
#'   for the interval-censored fit.
#'
#' @export
#' @rdname erf_methods
method_negbinomcens <- function(outcome = list(),
                                exposure = list(),
                                age = NULL,
                                population = NULL,
                                unit = NULL,
                                censor_lo = NULL,
                                censor_hi = NULL,
                                std_weights = NULL,
                                k_exposure = 20L,
                                k_conf = 10L,
                                em_iter = 8L,
                                ...) {
  rlang::check_installed("mgcv", "to use `method_negbinomcens()`.")

  arg::arg_list(outcome)
  arg::arg_list(exposure)

  arg::arg_count(k_exposure)
  arg::arg_gt(k_exposure, 1)
  arg::arg_count(k_conf)
  arg::arg_gt(k_conf, 1)
  arg::arg_count(em_iter)
  arg::arg_gt(em_iter, 0)

  out <- function(formula, data, weights = NULL, crossfit = 0L, n_grid = 100L,
                  cluster = NULL) {
    .c(Y, A, X) %=% process_formula(formula, data) #multi-assign

    n_obs <- length(A)

    if (is_null(weights)) {
      weights <- alloc(1.0, n_obs)
    }

    if (...length() > 0) {
      dots <- list(...)
      if ("cluster" %in% names(dots)) {
        stop("cluster is no longer an argument of method_negbinomcens(); supply it to estimate_erf().", call. = FALSE)
      }

      exposure <- utils::modifyList(exposure, dots)
      outcome <- utils::modifyList(outcome, dots)
    }

    if (is_null(exposure[["family"]])) {
      exposure$family <- gaussian()
    }

    # ---- Count-model columns (not part of the model formula) ----
    arg::when_not_null(age, arg::arg_name(data))
    arg::when_not_null(population, arg::arg_name(data))
    arg::when_not_null(unit, arg::arg_name(data))

    .age  <- if (is_null(age)) factor(alloc(1L, n_obs)) else factor(data[[age]])
    .pop  <- if (is_null(population)) alloc(1.0, n_obs) else as.numeric(data[[population]])
    .unit <- if (is_null(unit)) seq_len(n_obs) else data[[unit]]
    .clus <- if (is_null(cluster)) NULL else cluster[[1L]]

    # ---- Age-standardization weights ----
    alev <- levels(.age)
    if (is_null(std_weights)) {
      std_weights <-
        if (setequal(alev, c("65-74", "75-84", "85+")))
          c("65-74" = 18215 / 34834, "75-84" = 12360 / 34834, "85+" = 4259 / 34834)
      else setNames(rep(1 / length(alev), length(alev)), alev)
    }
    .wstd <- unname(std_weights[as.character(.age)])

    # ---- Interval-censoring bounds ----
    if (is_null(censor_lo) != is_null(censor_hi)) {
      arg::err("{.arg censor_lo} and {.arg censor_hi} must either both be specified or both be {.code NULL}")
    }

    arg::when_not_null(censor_lo, arg::arg_name(data))
    arg::when_not_null(censor_hi, arg::arg_name(data))

    if (is_null(censor_lo)) {
      if (anyNA(Y)) {
        arg::err("the outcome cannot contain {.val {NA}} values unless {.arg censor_lo} and {.arg censor_hi} give interval bounds")
      }
      .cens <- alloc(FALSE, n_obs)
      .cens_lo <- .cens_hi <- as.numeric(Y)
    }
    else {
      .cens_lo <- as.numeric(data[[censor_lo]])
      .cens_hi <- as.numeric(data[[censor_hi]])
      .cens <- is.na(Y) & is.finite(.cens_lo) & is.finite(.cens_hi)
    }

    # ---- Keep finite exposure, positive population, complete covariates, and an
    #      observed-or-censored outcome ----
    keep <- which(is.finite(A) & is.finite(.pop) & .pop > 0 &
                    stats::complete.cases(X) & (!is.na(Y) | .cens))

    Y <- as.numeric(Y)[keep]
    A <- A[keep]
    X <- ss(X, keep)
    .age <- droplevels(.age[keep])
    .pop <- .pop[keep]
    .unit <- .unit[keep]
    .wstd <- .wstd[keep]
    .cens <- .cens[keep]
    .cens_lo <- .cens_lo[keep]
    .cens_hi <- .cens_hi[keep]
    weights <- weights[keep]
    if (!is_null(.clus)) .clus <- .clus[keep]

    n <- length(A)

    # Working outcome; censored cells initialized to the interval midpoint.
    ywork <- Y
    ywork[.cens] <- (.cens_lo[.cens] + .cens_hi[.cens]) / 2

    # ---- Unit bookkeeping: units are the analysis rows for the ERF ----
    uf <- qF(.unit)
    ul <- levels(uf)
    n_unit <- length(ul)
    fr <- match(ul, as.character(.unit))          # a representative obs per unit

    A_unit <- A[fr]
    attr(A_unit, "varname") <- attr(A, "varname")
    A_vals <- seq(min(A_unit), max(A_unit), length.out = n_grid)

    # One cluster per unit. Taking the value at a representative observation is only
    # correct if the cluster does not vary within the unit, so check rather than assume:
    # a unit straddling two clusters would otherwise be silently assigned to one of them,
    # putting its cluster-mates in a different cross-fitting fold.
    if (!is_null(.clus)) {
      clus_per_unit <- fndistinct(as.character(.clus), g = uf)

      if (any(clus_per_unit > 1L)) {
        bad <- ul[clus_per_unit > 1L]
        arg::err("the clustering variable must not vary within a {.arg unit}, but it does for {cli::qty(length(bad))} unit{?s} ({.val {utils::head(bad, 5L)}})")
      }
    }

    clus_unit <- if (is_null(.clus)) ul else .clus[fr]

    # ---- Outcome GAM setup (age as factor; log(population) offset) ----
    cov_names <- names(X)
    is_cont <- vapply(cov_names, function(v) {
      is.numeric(X[[v]]) && fnunique(X[[v]]) > 10L
    }, logical(1L))
    
    # Get continuous and discrete variables
    cont_vars <- cov_names[is_cont]
    disc_vars <- cov_names[!is_cont]
    
    # Age term
    age_term <- if (!is_null(age) && nlevels(.age) > 1L) ".age" else NULL
    
    # Main effects: preserve the original specification.
    main_terms <- c(
      sprintf("s(A, bs = 'cr', k = %d)", k_exposure),
      if (length(cont_vars)) {
        sprintf("s(%s, bs = 'cr', k = %d)", cont_vars, k_conf)
      },
      disc_vars,
      age_term
    )
    
    # Smooth exposure-by-continuous-covariate interactions.
    interaction_terms <- c(
      if (length(cont_vars)) {
        sprintf(
          "ti(A, %s, bs = c('cr', 'cr'), k = c(%d, %d))",
          cont_vars,
          min(k_exposure, 5L),
          min(k_conf, 5L)
        )
      },
      
      # Parametric exposure interactions for discrete covariates.
      if (length(disc_vars)) {
        paste0("A:", disc_vars)
      },
      
      # Allow age-specific exposure slopes.
      if (length(age_term)) {
        "A:.age"
      }
    )
    
    # rhs <- paste(c(sprintf("s(A, bs = 'cr', k = %d)", k_exposure),
    #                if (any(is_cont)) sprintf("s(%s, bs = 'cr', k = %d)", cov_names[is_cont], k_conf),
    #                cov_names[!is_cont],
    #                age_term),
    #              collapse = " + ")
    
    rhs <- paste(c(main_terms, interaction_terms), collapse = " + ")
    
    out_form <- stats::as.formula(sprintf(".ywork ~ %s + offset(log(.pop))", rhs))

    dat_out <- cbind.data.frame(.ywork = ywork, A = A, .pop = .pop, .age = .age, X)

    # ---- GPS design (one row per unit) ----
    gps_data <- cbind.data.frame(A = A_unit, ss(X, fr))
    gps_form <- stats::reformulate(cov_names, response = "A")

    # ---- Helpers ----
    # E[y | lo <= y <= hi] under NB(mu, theta)
    cint <- function(mu, lo, hi, th) {
      vapply(seq_along(mu), function(i) {
        k <- lo[i]:hi[i]
        p <- stats::dnbinom(k, mu = mu[i], size = th)
        sum(k * p) / sum(p)
      }, numeric(1L))
    }

    # Censored-NB GAM via EM over mgcv::bam, fit on the training observations.
    fit_em <- function(id_train) {
      dtr <- ss(dat_out, id_train)
      w_tr <- weights[id_train]
      cen <- .cens[id_train]
      lo <- .cens_lo[id_train]
      hi <- .cens_hi[id_train]

      theta <- NULL
      m <- NULL
      th <- NA_real_

      for (it in seq_len(em_iter)) {
        args_bam <- utils::modifyList(
          outcome,
          list(formula = out_form, family = mgcv::nb(theta = theta),
               data = dtr, weights = w_tr, method = "fREML", discrete = TRUE))

        m <- do.call(mgcv::bam, args_bam)
        th <- m$family$getTheta(TRUE)

        if (!any(cen)) {
          break
        }

        mu <- as.numeric(stats::predict(m, type = "response"))
        enew <- cint(mu[cen], lo[cen], hi[cen], th)
        delta <- max(abs(enew - dtr$.ywork[cen]))
        dtr$.ywork[cen] <- enew

        if (delta < 1e-3) {
          break
        }

        theta <- th
      }

      attr(m, "theta") <- th
      m
    }

    # GPS conditional density (mirrors method_glm()'s dens_fun).
    make_dens_fun <- function(pimod) {
      fam <- family(pimod)$family
      if (fam == "Gamma") {
        alpha_hat <- 1 / summary(pimod)$dispersion
        function(a, fitted) dgamma(a, shape = alpha_hat, rate = alpha_hat / fitted)
      }
      else if (fam == "poisson") {
        function(a, fitted) dpois(a, fitted)
      }
      else {
        sd_hat <- sigma(pimod)
        function(a, fitted) dnorm(a, fitted, sd_hat)
      }
    }

    # ---- Fold assignment (by cluster when available, else by unit) ----
    if (crossfit <= 1L) {
      K <- 1L
      unit_fold <- alloc(1L, n_unit)
    }
    else {
      K <- as.integer(crossfit)
      unit_fold <- make_folds(K, n_unit, if (is_null(.clus)) NULL else clus_unit)
    }

    obs_unit <- match(as.character(.unit), ul)     # unit position for each obs

    # ---- Accumulate unit-level nuisance ----
    Y_unit <- numeric(n_unit)
    mu_hat <- numeric(n_unit)
    mu_hat_mat <- matrix(0.0, n_unit, n_grid)
    den_dens <- numeric(n_unit)
    pi_hat_mat <- matrix(0.0, n_unit, n_grid)
    thetas <- numeric(K)

    for (k in seq_len(K)) {
      eval_units <- if (K == 1L) seq_len(n_unit) else which(unit_fold == k)
      train_units <- if (K == 1L) seq_len(n_unit) else which(unit_fold != k)
      id_train <- if (K == 1L) seq_len(n) else which(obs_unit %in% train_units)
      id_eval <- if (K == 1L) seq_len(n) else which(obs_unit %in% eval_units)

      ## ---- Outcome model (observation level) ----
      m <- fit_em(id_train)
      thetas[k] <- attr(m, "theta")

      dat_e <- ss(dat_out, id_eval)
      pop_e <- .pop[id_eval]
      wstd_e <- .wstd[id_eval]
      uf_e <- qF(.unit[id_eval])

      mu_e <- as.numeric(stats::predict(m, newdata = dat_e, type = "response",
                                        discrete = FALSE))
      rate_e <- mu_e * 1e5 / pop_e

      yw_e <- Y[id_eval]
      ce <- .cens[id_eval]
      if (any(ce)) {
        yw_e[ce] <- cint(mu_e[ce], .cens_lo[id_eval][ce], .cens_hi[id_eval][ce],
                         attr(m, "theta"))
      }
      ywrate_e <- yw_e * 1e5 / pop_e

      mu_hat[eval_units] <- fmean(rate_e, g = uf_e, w = wstd_e)
      Y_unit[eval_units] <- fmean(ywrate_e, g = uf_e, w = wstd_e)

      ## ---- Outcome g-computation: mu_hat_mat[i, j] = E[rate | a_j, X_i] ----
      mm <- vapply(A_vals, function(a0) {
        rate_g <- as.numeric(stats::predict(m, newdata = transform(dat_e, A = a0),
                                            type = "response",
                                            discrete = FALSE)) * 1e5 / pop_e
        fmean(rate_g, g = uf_e, w = wstd_e)
      }, numeric(length(eval_units)))

      # # Log-additive shortcut: rate(a, X) = exp(s(a)) * B(X).
      # s_a <- as.numeric(stats::predict(m, newdata = dat_e, type = "terms", terms = "s(A)",
      #                                  discrete = FALSE))
      # s_a_unit <- fmean(s_a, g = uf_e)
      # Bunit <- mu_hat[eval_units] / exp(s_a_unit)
      # d1 <- ss(dat_e, 1L)
      # s_grid <- vapply(A_vals, function(a0) {
      #   d1$A <- a0
      #   as.numeric(stats::predict(m, newdata = d1, type = "terms", terms = "s(A)",
      #                             discrete = FALSE))
      # }, numeric(1L))
      # mm <- outer(as.numeric(Bunit), exp(s_grid))

      mu_hat_mat[eval_units, ] <- mm

      ## ---- GPS / exposure density (unit level) ----
      pimod <- do.call(stats::glm,
                       utils::modifyList(exposure,
                                         list(formula = gps_form,
                                              data = ss(gps_data, train_units))),
                       quote = TRUE)

      dfun <- make_dens_fun(pimod)
      fit_e <- stats::predict(pimod, newdata = ss(gps_data, eval_units), type = "response")

      den_dens[eval_units] <- dfun(A_unit[eval_units], fit_e)
      pi_hat_mat[eval_units, ] <- vapply(A_vals, function(a0) dfun(a0, fit_e),
                                         numeric(length(eval_units)))
    }

    nb_theta <- mean(thetas)

    pidens <- function(subset) {
      list(den_dens = ss(den_dens, subset),
           pi_hat_mat = ss(pi_hat_mat, subset))
    }

    mupred <- function(subset) {
      list(mu_hat = ss(mu_hat, subset),
           mu_hat_mat = ss(mu_hat_mat, subset))
    }

    cluster_df <- NULL
    if (!is_null(.clus)) {
      cluster_df <- data.frame(clus_unit) |>
        setNames(if (is_null(cluster)) "cluster" else names(cluster))
    }

    res <- list(Y = Y_unit,
                A = A_unit,
                A_vals = A_vals,
                pidens = pidens,
                mupred = mupred,
                weights = alloc(1.0, n_unit),
                cluster = cluster_df,
                # Unit-level covariate frame (one row per unit, in `A` order) so that
                # estimate_erf() can expose analysis-row-aligned modifier data for
                # downstream trees (e.g., honest_tree2(), rpart_erf()).
                data = ss(X, fr))

    attr(res, "nb_theta") <- nb_theta

    res
  }

  attr(out, "method_name") <- "negbinomcens"
  class(out) <- c(class(out), "erf_nuisance")

  out
}

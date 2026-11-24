#' Estimate an Exposure-Response Function
#'
#' @description
#' `estimate_erf()` estimates an exposure-response function (ERF) adjusting for measured covariates using the doubly robust estimator of Kennedy et al. (2017).
#'
#' @param formula a formula specifying the model. The should be in the form `Y ~ A | X1 + X2`, where `Y` is the outcome, `A` is the continuous exposure, and `X1` and `X2` are the covariates to adjust for.
#' @param data a data frame containing the variables in the model.
#' @param method an [`<erf_method>`][erf_methods] object specifying the method used to estimate the nuisance functions. The default is the output of [method_glm()], which fits generalized linear models for the exposure and outcome.
#' @param weights an optional numeric vector of sampling weights, with one value per row of `data`, or a string naming a numeric column of `data` to take them from. Default is `NULL` to weight every unit equally.
#' @param crossfit whether to crossfit the nuisance estimation. If `FALSE` or 1 (the default), the nuisance models will be trained and evaluated within the same sample. If an integer greater than 1, crossfitting will occur with the given number of folds. `TRUE` is the same as 2. When `cluster` is supplied, folds are formed over whole clusters. See Details.
#' @param cluster an optional one-sided formula (e.g., `~county`), vector, or one-column data frame giving the clustering variable, with one value per row of `data`. Only a single clustering variable is allowed. Clustering is used for two purposes: the cluster-robust variance of the ERF, and the cross-fitting folds, which are formed over whole clusters so that no observation's nuisance prediction is trained on its own cluster-mates. Default is `NULL` for no clustering.
#' @param range numeric; a numeric vector corresponding either to the lower and upper bounds of the exposure values for which to compute the ERF or a single number corresponding to the middle quantile of the exposure Default is 1 to use the full exposure range.
#' @param n integer specifying the number of equally spaced grid points on which to compute the ERF anchor points. Default is 51; higher numbers increase computation time and size of the resulting object but improve accuracy.
#' @param bw an optional bandwidth parameter. Default is `NULL`, which estimates the optimal bandwidth using leave-one-out cross-validation as described in Kennedy et al (2017). Other allowable arguments include a number giving the bandwidth or a function used to implement a rule of thumb, such as [bw.SJ()].
#' @param trim the quantile at which to trim the inverse-probability weights (the ratio of the marginal to the conditional exposure density) when forming the pseudo-outcomes. Where the conditional density is very small, this ratio can be enormous and a few units dominate the fit; capping it trades a little bias for a large variance reduction. Default is `0.99`, i.e., weights above their 99th percentile are capped. Set to `1` to disable trimming. Applies to all `method`s and is recomputed within each subgroup when growing a tree with [rpart_erf()].
#' @param n_grid an integer specifying the number of grid points to use to numerically evaluate integrals. Default is `100L`.
#' @param verbose A logical value indicating whether to print progress messages. Default is `FALSE`.
#' @param d the order of the derivative of the ERF to estimate. `0`, the default, estimates the ERF itself; `1` estimates its first derivative, i.e., the rate at which the expected outcome changes per unit of exposure at each anchor value. No other values are allowed. The degree of the local polynomial follows `d` as \eqn{p = d+1}, which keeps \eqn{p-d} odd as local polynomial theory recommends for estimating a \eqn{d}th derivative, and the bandwidth is selected for that same degree; an ERF and its derivative will therefore generally not share a bandwidth when `bw` is left `NULL`.
#' @param conley_se `logical`; whether to compute the covariance of the ERF as a Conley (1999) spatial HAC estimate, which allows the influence functions of nearby clusters to be correlated. The point estimates are unchanged. Requires `cluster`, and every cluster must sit at a single location. Default is `FALSE`.
#' @param coords when `conley_se = TRUE`, a character vector of length 2 naming the variables in `data` holding longitude and latitude, in that order and in decimal degrees. Default is `c("longitude", "latitude")`.
#' @param cutoff when `conley_se = TRUE`, the distance in kilometers beyond which two clusters are treated as independent. A pair of clusters a distance \eqn{d} apart is weighted by the Bartlett kernel \eqn{1 - d/\mathrm{cutoff}}{1 - d/cutoff}, so the variance is understated unless `cutoff` is several times the distance over which the dependence persists. Default is 500.
#' @param ... Additional arguments passed to the underlying estimation function.
#'
#' @return
#' An `<erf>` object, which contains the following components:
#' \item{`erf`}{a function with two arguments:
#' \describe{
#'   \item{`value`}{the values of the exposure over which to evaluate the ERF}
#'   \item{`vcov`}{`logical`; whether to compute the covariance matrix at each value (default `TRUE`)}
#' }
#' This function returns a list of two elements:
#' \describe{
#'   \item{`est`}{the ERF estimate at the points supplied to `values`}
#'   \item{`vcov`}{when `vcov = TRUE`, the covariance matrix of the estimates}
#' }
#' }
#' \item{`formula`}{the model formula}
#' \item{`data`}{the original dataset used to estimate the ERF, including only the columns used in the estimation}
#' \item{`treat`}{the exposure variable, which has attribute `"varname"` with its name}
#' \item{`method`}{the `<erf_method>` object used to estimate the nuisance functionals}
#' \item{`fit`}{a list containing the output of `method` applied to the dataset}
#' \item{`weights`}{a vector of sampling weights; set to `1` if `NULL` in the original call to `erf()`}
#' \item{`bw`}{the bandwidth the ERF was smoothed with; always a number, whether `bw` was supplied as one, as a function, or left `NULL` to be selected by cross-validation}
#' \item{`values`}{a vector of anchor points determined by `range` and `n`}
#' \item{`trim`}{the quantile the density ratio was trimmed at}
#' \item{`d`}{the order of the derivative the ERF estimates; `0` for the ERF itself}
#' \item{`cluster`}{the clustering variable aligned to the analysis rows, or `NULL` when unclustered}
#' \item{`call`}{the original function call}
#'
#' @details
#' `estimate_erf()` estimates the ERF using the methods described in Kennedy et al. (2017). This involves estimating the nuisance functionals used to compute the pseudo-outcomes and fitting a local linear regression of the pseudo-outcomes on the exposure. The nuisance functionals include estimates of the conditional density of exposure and estimates of the outcome mean, estimated using the method supplied to `method`. The covariance matrix for the ERF estimates is computed using the efficient influence function.
#'
#' ## Crossfitting
#'
#' Crossfitting enhances the validity of the inferences made using the efficient influence function. First, the data is split into \eqn{k} independent folds, where \eqn{k} corresponds to the value supplied to `crossfit`. One fold is designated the evaluation fold, and the other \eqn{k-1} folds are combined and then used to fit the nuisance models. The nuisance estimates are evaluated in the evaluation fold. The process repeats, with each of the \eqn{k} folds eventually used as the evaluation fold. The nuisance estimates are then used to compute the pseudo-outcomes.
#'
#' @seealso [plot.erf()]; [rpart_erf()]
#'
#' @references
#' Kennedy, E. H., Ma, Z., McHugh, M. D., & Small, D. S. (2017). Non-parametric methods for doubly robust estimation of continuous treatment effects. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 79(4), 1229–1245. \doi{10.1111/rssb.12212}
#'
#' Conley, T. G. (1999). GMM estimation with cross sectional dependence. *Journal of Econometrics*, 92(1), 1–45. \doi{10.1016/S0304-4076(98)00084-0}
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
#' # Exposure on the left of the `|`, confounders on the right
#' e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 25L)
#' e
#'
#' # The ERF at chosen exposure values, with a covariance matrix
#' est <- e$erf(c(0, 1, 2), vcov = TRUE)
#' est$est
#' sqrt(diag(est$vcov))
#'
#' \donttest{
#' plot(e)
#'
#' # Cross-fit the nuisance models over 5 folds
#' e2 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 25L, crossfit = 5)
#'
#' # Clustered data: `cluster` gives both cluster-robust standard errors and
#' # cross-fitting folds formed over whole clusters
#' d$id <- factor(rep(seq_len(50), each = 10))
#' e3 <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 25L, cluster = ~id, crossfit = 5)
#' }
#'
#' @export
estimate_erf <- function(formula, data, method = method_glm(), weights = NULL, bw = NULL,
                         crossfit = FALSE, cluster = NULL, n = 51L, range = 1, n_grid = 100L,
                         trim = 0.99, verbose = FALSE, d = 0, 
                         conley_se = FALSE, 
                         coords = c("longitude", "latitude"), 
                         cutoff = 500, ...) {
  call <- match.call()

  # Process inputs
  A <- process_exposure(formula, data)

  # Process method; should be a function
  method <- process_nuissance(method)

  # Process weights
  weights <- process_weights(weights, data)

  # Process crossfit
  crossfit <- process_crossfit(crossfit)

  # Process trim
  trim <- process_trim(trim)

  # Process d
  d <- process_d(d)

  # Process cluster. Resolved before the nuisance models are fit, because clustering
  # determines the cross-fitting folds as well as the standard errors: whole clusters are
  # kept together in a fold so that no unit's nuisance prediction is trained on its own
  # cluster-mates.
  cluster <- resolve_cluster(cluster, data)
  
  # Process coords before the nuisance models are fit, so a bad one fails immediately.
  arg::arg_flag(conley_se)
  xy <- if (conley_se) process_coords(coords, cutoff, data, cluster)
  
  # Call method on inputs to estimate nuisances. A method may aggregate the data
  # (e.g., method_negbinomcens() collapses observations to units), in which case the
  # ERF, its bandwidth, weights, and clusters all live at the aggregated level. The
  # exposure, weights, and clusters used below are therefore taken from `fit`.
  fit <- fit_nuisance(formula, data, method, weights, crossfit, n_grid,
                      cluster = cluster, verbose = verbose)
  
  if (!is_null(xy)) {
    fit[["coords"]] <- xy[match(as.character(fit[["cluster"]][[1L]]),
                                as.character(cluster[[1L]])), , drop = FALSE]
    fit[["cutoff"]] <- cutoff
  }
  
  A_erf <- fit[["A"]]
  if (is_null(attr(A_erf, "varname"))) {
    attr(A_erf, "varname") <- attr(A, "varname")
  }
  weights_erf <- fit[["weights"]]

  # Process bw
  bw <- process_bw(bw, A = A_erf)

  values <- process_range(range, n, A = A_erf, w = weights_erf)

  # Every method returns the clustering aligned to its own analysis rows -- units when
  # the method aggregates, observations otherwise -- so the variance clustering comes
  # from `fit` unconditionally.
  cluster_erf <- process_cluster(cluster = fit[["cluster"]],
                                 data = fit[["cluster"]])

  # Pass to internal function
  erf <- estimate_erf_internal(fit = fit, values = values, vcov = TRUE,
                               bw = bw, cluster = cluster_erf, verbose = verbose,
                               d = d, trim = trim)

  if (verbose) {
    message("Done!")
  }

  # `data` is stored aligned to the analysis rows (i.e. to `fit$A`). A method that
  # aggregates observations to units (e.g. method_negbinomcens()) returns its
  # unit-level covariate frame in `fit$data`; otherwise the observation-level model
  # frame is used, unchanged.
  out <- list(erf = erf,
              formula = formula,
              data = fit[["data"]] %or% data[names(data) %in% get_varnames(formula)],
              treat = A_erf,
              method = method,
              fit = fit,
              weights = weights_erf,
              bw = attr(erf, "bw"),
              values = values,
              trim = trim,
              # Stored so that update(bw = ) can re-smooth at the same coefficient order
              # rather than silently falling back to the ERF itself.
              d = d,
              # Aligned to the analysis rows, i.e. to `treat`. NULL when unclustered.
              cluster = fit[["cluster"]],
              call = call)

  class(out) <- "erf"

  out
}

estimate_erf_internal <- function(fit, values, subset = NULL, bw = NULL,
                                  vcov = TRUE, cluster = NULL, verbose = FALSE, d = 0,
                                  trim = 1) {

  # Estimate Pseudo-Y and related components from nuisance
  res <- estimate_pseudo_Y(fit, subset = subset, trim = trim)

  # The local polynomial order follows the coefficient being returned: p = d + 1 keeps
  # p - d odd, the configuration local polynomial theory recommends for estimating the
  # d-th derivative (Fan & Gijbels, 1996). The bandwidth is selected for the same order
  # that is then fit.
  p <- d + 1L
  
  coords <- fit[["coords"]]
  if (!is_null(coords) && !is_null(subset)) {
    coords <- coords[subset, , drop = FALSE]
  }
  
  # Compute bw. `res` is already restricted to `subset`, so the bandwidth is selected on
  # the exposures and weights actually being smoothed -- not on the full sample, which
  # would not even be the same length as `pseudo_Y`.
  bw <- get_bw(A = res[["A"]], Y = res[["pseudo_Y"]],
               weights = res[["weights"]], bw = bw, p = p,
               verbose = verbose)

  # Compute ERF from pseudo outcomes
  make_erf(res, anchor_values = values, bw = bw, vcov = vcov,
           cluster = cluster, p = p, verbose = verbose, d = d,
           coords = coords, cutoff = fit[["cutoff"]])
}

fit_nuisance <- function(formula, data, method, weights, crossfit, n_grid, cluster = NULL,
                         verbose = TRUE) {
  if (verbose) {
    message("Fitting nuisance models...")
  }

  method(formula = formula, data = data,
         weights = weights, crossfit = crossfit,
         n_grid = n_grid, cluster = cluster)
}

estimate_pseudo_Y <- function(fit, subset = NULL, trim = 1) {

  Y <- fit[["Y"]]
  A <- fit[["A"]]
  A_vals <- fit[["A_vals"]]
  weights <- fit[["weights"]]

  n <- length(A)

  if (is_null(subset)) {
    subset <- rep.int(TRUE, n)
  }

  pred_ID <- seq_len(n)[subset]

  if (is_null(weights)) {
    weights <- rep.int(1.0, n)
  }

  weights <- weights[pred_ID]

  ## Evaluated conditional density function
  pidens <- fit[["pidens"]](pred_ID)

  ### Denominator density
  den_dens <- pidens$den_dens

  ### Matrix of conditional density for each X * A_grid grid unit
  pi_hat_mat <- pidens$pi_hat_mat

  ### Marginal density:
  num_dens <- predict(smooth.spline(A_vals, fmean(pi_hat_mat, w = weights)),
                      x = A[pred_ID])$y


  ## Outcome prediction function
  mupred <- fit[["mupred"]](pred_ID)

  ## Predicted Y for each unit
  mu_hat <- mupred$mu_hat

  ## Matrix of predicted Y for each X * A_grid grid unit
  mu_hat_mat <- mupred$mu_hat_mat

  ## Integral of predicted Y over X for each unit's value of A
  m_hat <- predict(smooth.spline(A_vals, fmean(mu_hat_mat, w = weights)),
                   x = A[pred_ID])$y

  ## Inverse-probability weights: marginal over conditional exposure density
  ipw <- num_dens / den_dens

  ## Optionally trim the weights to stabilize the pseudo-outcomes. Where the
  ## conditional density is very small (sparse regions of the exposure
  ## distribution), the ratio can be enormous and a handful of units then dominate
  ## the local-linear fit. Trimming caps the ratio at its `trim` quantile, trading a
  ## little bias for a large variance reduction. This is applied here, rather than in
  ## the nuisance methods, so that it (a) applies to every method, (b) bounds the
  ## ratio that is actually used, and (c) is recomputed within each `subset`, so
  ## subgroup ERFs are stabilized too. `num_dens` comes from a smoothing spline and
  ## can dip below 0 in sparse tails; negative weights would flip the sign of the
  ## residual correction, so they are floored at 0.
  ipw <- pmax(ipw, 0)

  if (trim < 1) {
    finite_ipw <- ipw[is.finite(ipw)]

    if (!is_null(finite_ipw)) {
      cap <- unname(.quantile(finite_ipw, w = weights[is.finite(ipw)], probs = trim))

      ipw[!is.finite(ipw)] <- cap
      ipw <- pmin(ipw, cap)
    }
  }
  else {
    ipw[!is.finite(ipw)] <- 0
  }

  ## Pseudo outcome for each unit
  pseudo_Y <- (Y[pred_ID] - mu_hat) * ipw + m_hat

  list(pseudo_Y = pseudo_Y,
       mu_hat_mat = mu_hat_mat,
       pi_hat_mat = pi_hat_mat,
       ipw = ipw,
       A = A[pred_ID],
       # The exposure grid the columns of `pi_hat_mat` and `mu_hat_mat` are indexed by.
       # Retained so that consumers averaging those matrices over units (e.g. the "gps"
       # integral weights in get_int_weights()) can interpolate on the correct abscissa.
       A_vals = A_vals,
       weights = weights)
}

# The function an <erf> exposes as its `erf` component: it maps requested exposure values
# onto the estimates (and covariance) held at the anchor values, interpolating between
# anchors where necessary. Factored out of make_erf() so that the thin per-node ERFs
# stored on an <erf_tree> interpolate in exactly the same way.
erf_interpolator <- function(anchor_values, est, vcov_mat = NULL) {
  force(anchor_values)
  force(est)
  force(vcov_mat)

  function(values, vcov = TRUE) {
    val_mat <- matrix(0, nrow = length(values),
                      ncol = length(anchor_values))

    matches <- max.col(-abs(outer(values, anchor_values, "-")), "first")

    to_interp <- !check_if_zero(values - anchor_values[matches])

    if (any(to_interp)) {
      val_mat[to_interp, ] <- get_locpoly_w(x = values[to_interp],
                                            v = anchor_values)
    }

    if (!all(to_interp)) {
      val_mat[cbind(which(!to_interp), matches[!to_interp])] <- 1
    }

    list(
      exposure = values,
      est = drop(val_mat %*% est),
      vcov = if (vcov && !is_null(vcov_mat)) val_mat %*% vcov_mat %*% t(val_mat)
    )
  }
}

# Returns a function that produces the ERF (which is a function)
make_erf <- function(res, anchor_values, bw = NULL, vcov = TRUE, cluster = NULL,
                     p = d + 1L, verbose = TRUE, d = 0,
                     coords = NULL, cutoff = NULL) {
  if (verbose) {
    message("Computing ERF...")
  }

  k <- dnorm

  # The design matrix is built on the bandwidth-standardized exposure (A - a)/bw, so the
  # fitted coefficient of order d estimates theta^(d)(a) * bw^d / d! rather than the
  # derivative itself. Undoing that here puts both the estimate and, below, the influence
  # function on the scale of theta^(d)(a). For d = 0 the factor is exactly 1.
  d_scale <- factorial(d) / bw^d

  A0 <- res[["A"]]
  pseudo_Y0 <- res[["pseudo_Y"]]
  weights <- res[["weights"]]

  n_grid <- ncol(res[["mu_hat_mat"]])

  est <- numeric(length(anchor_values))

  if (vcov) {
    n <- length(A0)

    m_hat <- fmean(res[["mu_hat_mat"]], w = weights)

    p_hat <- fmean(res[["pi_hat_mat"]], w = weights)

    Imat <- (res[["mu_hat_mat"]] %r-% m_hat) %r*% p_hat

    A_vals <- seq(min(A0), max(A0), length.out = n_grid) #density grid

    inf <- vector("list", length(anchor_values))
  }

  for (i in seq_along(anchor_values)) {
    A_std <- (A0 - anchor_values[i]) / bw

    G <- get_G(A_std, p)
    K <- k(A_std) / bw

    # ERF is the intercept for each local poly reg
    fiti <- lm.wfit(x = G, y = pseudo_Y0, w = K * weights)

    est[i] <- fiti$coefficients[d + 1L] * d_scale

    if (vcov) {
      A_grid_std <- (A_vals - anchor_values[i]) / bw
      G_grid <- get_G(A_grid_std, p)
      K_grid <- k(A_grid_std) / bw

      r <- fiti$residuals

      V <- do.call("rbind", lapply(seq_len(ncol(G)), function(j) {
        intfni.mat <- Imat %r*% (G_grid[, j] * K_grid)

        # Trapezoidal approximation to integral
        int <- rowSums(((intfni.mat[, -1L, drop = FALSE] + intfni.mat[, -n_grid, drop = FALSE]) / 2) %r*% diff1(A_vals))

        weights * (G[, j] * K * r + int)
      }))

      ## Influence function

      ### Compute the weighted design matrix cross-product
      des_matrix <- crossprod(G, G * (K * weights / n))

      ### Get the robust pseudo-inverse matrix
      p_inv <- pseudo_inverse(des_matrix)

      ### Multiply by V and extract row (d + 1)
      inf_matrix <- p_inv %*% V
      inf[[i]] <- inf_matrix[d + 1L, ] * d_scale / sqrt(n)
    }
  }

  anchors <- list(est = est)

  if (vcov) {
    inf <- do.call("cbind", inf)
    
    if (!is_null(coords)) {
      anchors$vcov <- conley_vcov(inf, coords[, 1L], coords[, 2L], cutoff)
    }

    else if (is_null(cluster)) {
      anchors$vcov <- cov(inf)
    }
    
    else {
      # The influence rows in `inf` are scaled by 1/sqrt(n) so that the iid path
      # `cov(inf)` (which divides by n - 1) is correctly calibrated. The clustered
      # sum of squares must therefore also be divided by n to be on the same scale;
      # with the small-sample factor `adj` = G / (G - 1), a singleton-cluster
      # clustering then reduces to `cov(inf)`.
      anchors$vcov <- Reduce("+", lapply(seq_along(cluster[["cluster"]]), function(j) {
        crossprod(fsum(inf, g = cluster[["cluster"]][[j]])) * cluster[["adj"]][j]
      })) / n
    }
  }

  erf_fun <- erf_interpolator(anchor_values, anchors$est, anchors$vcov)

  attr(erf_fun, "bw") <- bw
  attr(erf_fun, "res") <- res

  erf_fun
}

# Returns bandwidth (number)
get_bw <- function(A, Y, weights, bw = NULL, p = 1, verbose = TRUE) {
  if (!is_null(bw) && is.numeric(bw)) {
    return(bw)
  }

  n <- length(A)

  if (is_null(bw)) {
    # leave-one-out cross-validation to select bandwidth
    if (verbose) {
      message("Estimating optimal bandwidth...")
    }

    k <- dnorm

    LOOCV_obj <- get_LOOCV_obj(p, A, Y, k, weights)

    opt <- optimize(LOOCV_obj, interval = c(1 / length(A), 3) * sd(A))

    bw <- opt$minimum
  }
  else if (is.function(bw)) {
    if (all(c("x", "y") %in% rlang::fn_fmls_names(bw))) {
      bw <- bw(x = A, y = Y)
    }
    else {
      bw <- bw(A)
    }
  }

  arg::arg_number(bw)
  arg::arg_gt(bw, 0)

  bw
}

# Construct polynomial basis of `a`
get_G <- function(a, p = 1L) {
  if (p == 0L) {
    return(matrix(1, nrow = length(a)))
  }

  if (p == 1L) {
    return(cbind(1, a))
  }

  cbind(1, outer(a, seq_len(p), "^"))
}

get_locpoly_w <- function(x, v, p = 3L) {
  n <- length(v)

  bw <- mean(diff1(v)) / 2

  w <- dnorm(outer(x, v, "-"), sd = bw)

  rs <- rowSums(w)

  zeros <- check_if_zero(rs)

  if (any(zeros)) {
    w[zeros, ] <- 0
  }

  if (!all(zeros)) {
    w[!zeros, ] <- w[!zeros, ] / rs[!zeros]
  }

  if (p == 0) {
    return(w)
  }

  ind <- seq_len(n)

  v0 <- c(v, x)

  dv0 <- cbind(1, poly(v0, degree = p, raw = FALSE, simple = TRUE))

  dx <- ss(dv0, -ind)

  lw <- matrix(0, nrow = length(x), ncol = n)

  for (i in seq_along(x)) {
    wi <- w[i, ]
    pos <- which(wi > 1e-10)

    if (is_null(pos)) {
      arg::err("all kernel weights estimated as 0, indicating severe extrapolation")
    }

    wi <- wi[pos]
    wi <- wi / sum(wi)

    lw[i, pos] <- drop(ss(dx, i) %*% solve(crossprod(ss(dv0, pos) * sqrt(wi)), t(ss(dv0, pos) * wi)))
  }

  lw
}

get_LOOCV_obj <- function(p = 1L, A, Y, k, weights, max_size = 5e3) {
  arg::arg_count(max_size)

  if (length(A) > max_size) {
    sub <- sample.int(length(A), max_size)
  }
  else {
    sub <- seq_along(A)
  }

  if (p > 1L) {
    function(h) {
      res_h <- vapply(sub, function(i) {
        a0_std <- (A[i] - A[-i]) / h

        G <- get_G(a0_std, p)

        wsq <- sqrt(weights[-i] * k(a0_std) / h)

        Y[i] - .lm.fit(wsq * G, wsq * Y[-i])$coefficients[1L]
      }, numeric(1L))

      mean(res_h^2)
    }
  }
  else if (p == 1L) {
    Dm  <- -outer(A[sub], A, "-")

    function(h) {
      W   <- k(Dm / h) / h

      s0  <- rowSums(W)
      s1  <- rowSums(W * Dm)
      s2  <- rowSums(W * Dm^2)
      Wdy <- rowSums((W * Dm) %r*% Y)

      denom <- s0 * s2 - s1^2
      yhat  <- (s2 * (W %*% Y) - s1 * Wdy) / denom
      S_ii  <- (k(0) / h) * (s2 / denom)

      mean(((Y[sub] - yhat) / (1 - S_ii))^2)
    }
  }
  else {
    D  <- outer(A[sub], A, "-")

    function(h) {
      W <- k(D / h) / h

      s0  <- rowSums(W)

      yhat <- (W %*% Y) / s0
      S_ii <- k(0) / (h * s0)

      mean(((Y[sub] - yhat) / (1 - S_ii))^2)
    }
  }
}

# Conley covariance across anchor values, with Bartlett weights on great-circle distance
# in km. `inf` is the influence matrix of make_erf(), rows scaled by 1 / sqrt(n), and the
# scaling matches its clustered branch, which a vanishing cutoff reproduces.
conley_vcov <- function(inf, lon, lat, cutoff) {
  loc <- paste(lon, lat)
  id <- match(loc, unique(loc))
  
  psi <- rowsum(inf, id, reorder = FALSE)
  first <- !duplicated(id)
  
  phi <- lat[first] * pi / 180
  lambda <- lon[first] * pi / 180
  
  G <- nrow(psi)
  
  a <- sin(outer(phi, phi, "-") / 2) ^ 2 +
    outer(cos(phi), cos(phi)) *
    sin(outer(lambda, lambda, "-") / 2) ^ 2
  
  dist_km <- 2 * 6371 * asin(sqrt(pmin(a, 1)))
  W <- (dist_km < cutoff) * (1 - dist_km / cutoff)
  V <- crossprod(psi, W %*% psi) * (G / (G - 1)) / nrow(inf)
  V <- (V + t(V)) / 2
  
  if (!all(is.finite(V))) {
    return(V)
  }
  
  e <- eigen(V, symmetric = TRUE)
  pos <- e$values > 0
  
  if (all(pos)) {
    return(V)
  }
  
  Q <- e$vectors[, pos, drop = FALSE]
  
  Q %*% (e$values[pos] * t(Q))
}

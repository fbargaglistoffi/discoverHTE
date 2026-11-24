#' Recursive Partitioning on an Exposure-Response Function
#'
#' @description
#' Partitions a sample into subgroups by maximizing the difference between subgroup ERFs.
#'
#' @param erf an `<erf>` object; the output of a call to [estimate_erf]. When it was fit with `d = 1`, the subgroup ERFs are estimated at that same order, so the tree searches for differences in the derivative of the ERF rather than in the ERF itself.
#' @param modifiers a one-sided formula containing the potential effect modifiers, e.g., `~ M1 + M2`. If `NULL`, the default, will use all covariates present in `erf`. A unit missing the modifier a split is made on is sent to the right branch of that split, so the candidate cut points are read off the observed values alone.
#' @param metric a string corresponding to the name of the objective function to maximize at each split. Allowable options include `"L2"` for the \eqn{L_2} metric, `"L1"` for the \eqn{L_1} metric, `"fANOVA"` for the functional ANOVA metric, and `"AUC"` for the area-under-curve (AUC) metric. See Details for more information. Default is `"L2"`.
#' @param min_samples the smallest size for each leaf node. Default is 30.
#' @param max_depth the maximum depth of the trees. Default is 3.
#' @param int_weight_method a string corresponding to the name of the method used to weight the objective function. Allowable options include `"density"` to weight by the exposure density, `"gps"` to weight by the conditional exposure density, and `"none"` to weight each area of the ERF uniformly. See Details for more information. Default is `"density"`.
#' @param cont_var_splits how many levels continuous modifiers should be split into in each partition. Default is 5.
#' @param split_type string; how continuous modifiers should be split. Allowable options include `"quantile"`, which splits continuous modifiers into `cont_var_splits` quantiles (i.e., containing approximately the same number of units in each), and `"range"`, which splits continuous modifiers into evenly spaced bins on the scale of the modifier. Default is `"quantile"`.
#' @param erf_est an optional ERF fitted on an independent estimation sample, using the same formula and effect modifiers. When supplied, the splits learned on `erf` are applied to this sample to estimate the node curves.
#'
#' @return
#' An `<erf_tree>` object, which contains the following components:
#' \item{`tree`}{the recursive partitioning tree; see below}
#' \item{`M`}{a data frame containing the potential modifiers searched}
#' \item{`values`}{the anchor values of the ERF}
#' \item{`modifiers`}{the formula supplied to `modifiers`, or `NULL` when all covariates were searched}
#' \item{`d`}{the order of the derivative the ERF estimates; `0` for the ERF itself}
#' \item{`info`}{a list of the settings the tree was grown with: `metric`, `int_weight_method`, `max_depth`, and `min_samples`}
#'
#' The `tree` component contains the recursive partitioning tree. Each leaf (i.e., terminal) node is a list with the following components:
#'
#' \item{`type`}{"leaf" for leaf nodes}
#' \item{`path`}{a character vector identifying to which split of the modifiers the node belongs}
#' \item{`size`}{the number of units in the node}
#' \item{`reason`}{the reason for termination at the node}
#' \item{`erf`}{the ERF function for that node; see the `erf` component of the output of [estimate_erf()]}
#'
#' Each internal node is a list with the following components:
#' \item{`type`}{"internal" for internal nodes}
#' \item{`path`}{a character vector identifying to which split of the modifiers the node belongs}
#' \item{`size`}{the number of units in the node}
#' \item{`split_var`}{the name of the variable split at the node}
#' \item{`split_levels`}{the levels of the variable split at the node}
#' \item{`left`, `right`}{each is a list containing a leaf or internal node corresponding to the split occurring at the node}
#' \item{`metric`}{the value of the objective function for the split of the node}
#' \item{`erf`}{the ERF function for that node; see the `erf` component of the output of [estimate_erf()]}
#'
#' @details
#' `rpart_erf()` recursively splits the sample into subgroups that maximize a measure of effect heterogeneity in the ERFs. The measure is determined by the argument to `metric`. Below are the options it can take and the formulas used to compute them, where \eqn{\mathcal{A}} is the range of exposure values, \eqn{\hat\theta_L(a)} and \eqn{\hat\theta_R(a)} are the left and right branch ERFs, \eqn{\hat\theta(a)} is the parent ERF, and \eqn{w(a)} is the integral weighting function:
#' | `metric` | Formula |
#' |:---------|:--------|
#' |`"L2"`    | \eqn{\int_{\mathcal{A}}\left(\hat\theta_L(a) - \hat\theta_R(a)\right)^2 w(a) da} |
#' |`"L1"`    | \eqn{\int_{\mathcal{A}}\left| \hat\theta_L(a) - \hat\theta_R(a) \right| w(a) da}
#' |`"fANOVA"`| \eqn{p\int_{\mathcal{A}}\left(\hat\theta_L(a) - \hat\theta(a)\right)^2 w(a) da + (1-p)\int_{\mathcal{A}}\left(\hat\theta_R(a) - \hat\theta(a)\right)^2 w(a) da} |
#' |`"AUC"`   | \eqn{\left| \int_{\mathcal{A}} \hat\theta_L(a) w(a) da - \int_{\mathcal{A}} \hat\theta_R(a) w(a) da \right|}|
#'
#' In the `"fANOVA"` formula, \eqn{p} is the left branch's share of the node, \eqn{n_L/(n_L+n_R)}, computed from the weighted counts. Weighting each child's deviation by its share keeps a small subgroup that happens to deviate wildly from scoring like a large one, and makes `"fANOVA"` a balance-penalized version of `"L2"`: at the population level the two differ by a factor of \eqn{p(1-p)}, which is largest at an even split.
#'
#' The integral weighting function \eqn{w(a)} is determined by the argument to `int_weight_method`. In each case \eqn{\hat\pi_L(a)} and \eqn{\hat\pi_R(a)} are estimates of the exposure density in the left and right branches, and the weight is normalized to integrate to 1 over \eqn{\mathcal{A}}:
#' | `int_weight_method` | Formula |
#' |:---------|:--------|
#' |`"none"`    | \eqn{w(a) \propto 1} |
#' |`"density"`    | \eqn{w(a) \propto \hat\pi_L(a)\hat\pi_R(a)}, with each \eqn{\hat\pi} the weighted marginal (kernel) density of the exposure in that branch |
#' |`"gps"`| \eqn{w(a) \propto \hat\pi_L(a)\hat\pi_R(a)}, with each \eqn{\hat\pi} the branch-averaged *conditional* exposure density \eqn{n_L^{-1}\sum_{i \in L}\hat\pi(a \mid X_i)} |
#'
#' Both `"density"` and `"gps"` therefore weight by the product of the two branches' exposure densities, which concentrates the objective where *both* subgroups have exposure mass and so keeps the metric from being driven by extrapolation into a region only one subgroup occupies. Because \eqn{w(a)} depends on the two branches, it is recomputed for each candidate split. Where neither branch has support the weight is 0; if no exposure value has mass in both branches, the anchor values are weighted uniformly instead.
#'
#' The ERF at each node is estimated from the nuisance functionals stored in `erf`, with the pseudo-outcomes recomputed within the node, at the bandwidth of `erf` and with its trimming and its clustering, so the node bands are cluster-robust when the ERF is. They are estimated on the same sample the tree was grown on, however, and so do not account for the search over candidate splits that selected the nodes; [honest_erf_tree()] grows the tree on one part of the sample and estimates the node ERFs on the other for that reason.
#'
#' @seealso [estimate_erf()]; [honest_erf_tree()] for a tree whose node ERFs are estimated on an independent sample; [summary.erf_tree()], [plot.erf_tree()], and [`[[.erf_tree`] for reading the tree
#'
#' @references
#' Kennedy, E. H., Ma, Z., McHugh, M. D., & Small, D. S. (2017). Non-parametric methods for doubly robust estimation of continuous treatment effects. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 79(4), 1229–1245. \doi{10.1111/rssb.12212}
#'
#'
#' @examples
#' # Simulated data with an effect modifier: the exposure-response slope
#' # differs by the binary modifier M.
#' set.seed(1234)
#' n <- 500
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' M <- factor(rbinom(n, 1, 0.5))
#' A <- 1 + 0.5 * X1 + rnorm(n)
#' Y <- A * (1 + 1.5 * (M == "1")) + X1 + X2 + rnorm(n)
#' d <- data.frame(Y, A, X1, X2, M)
#'
#' e <- estimate_erf(Y ~ A | X1 + X2 + M, data = d, n = 25L)
#'
#' # Search M and X1 for effect modification
#' tr <- rpart_erf(e, ~ M + X1, max_depth = 2, min_samples = 50)
#' tr
#' summary(tr)
#'
#' # The ERF at a single node, by its position in summary()
#' plot(tr[[2]])
#'
#' \donttest{
#' # The whole tree, each internal node showing its two child ERFs
#' plot(tr)
#' }
#'
#' @export
rpart_erf <- function(erf, modifiers = NULL, metric = "L2",
                      min_samples = 30, max_depth = 3, int_weight_method = "density",
                      cont_var_splits = 5L, split_type = "quantile",
                      erf_est = NULL) {

  # Check erf. Note a node ERF inherits from "erf" but carries neither the data nor the
  # fitted nuisance models that growing a tree requires.
  arg::arg_is(erf, "erf")
  check_full_erf(erf)

  # Process modifiers; formula input, model.frame output
  M <- process_modifiers(modifiers, erf)

  # Check int_weight_method
  int_weight_method <- process_int_weight_method(int_weight_method)

  # Process objective function
  metric <- process_metric(metric)

  # Check min_samples
  arg::arg_count(min_samples)
  arg::arg_gt(min_samples, 1)

  # Check max_depth
  arg::arg_count(max_depth)

  arg::arg_count(cont_var_splits)
  split_type <- arg::match_arg(split_type, c("quantile", "range"))

  # Build tree
  tree <- build_tree(M = M,
                     fit = erf$fit,
                     bw = erf$bw,
                     values = erf$values,
                     treat = erf$treat,
                     max_depth = max_depth,
                     min_samples = min_samples,
                     metric = metric,
                     int_weight_method = int_weight_method,
                     cont_var_splits = cont_var_splits,
                     split_type = split_type,
                     erf_parent = erf$erf,
                     trim = erf$trim %or% 1,
                     # The subgroup ERFs are estimated at the same coefficient order as
                     # the ERF being partitioned, so a tree grown on a derivative ERF
                     # splits on differences in the derivative.
                     d = erf$d %or% 0,
                     # The node ERFs carry the same clustering as the ERF being
                     # partitioned, so their bands are cluster-robust when it is.
                     cluster = erf$cluster)

  out <- list(tree = tree,
              M = M,
              method = erf$method,
              values = erf$values,
              treat = erf$treat,
              treat_name = attr(erf$treat, "varname"),
              pseudo_Y = attr(erf$erf, "res")[["pseudo_Y"]],
              weights = erf$weights,
              # The modifiers formula and the coefficient order are kept so that
              # refit_erf_tree() can resolve the same modifiers against another ERF and
              # check that it was estimated at the same order.
              modifiers = modifiers,
              d = erf$d %or% 0,
              info = list(
                metric = metric,
                int_weight_method = int_weight_method,
                max_depth = max_depth,
                min_samples = min_samples
              )
  )

  class(out) <- "erf_tree"
  
  # Add erf_est, which is fitted independently from erf, to allow for honesty.
  # If not NULL, rules are first discovered with erf, then applied on erf_est.
  # erf and erf_est MUST share the same formula and effect modifiers.
  if (!is.null(erf_est)) {
    out <- apply_tree(out, erf_est)
  }

  out
}

build_tree <- function(M, fit, split = rep.int(TRUE, nrow(M)), depth = 0, path = "root",
                       bw = NULL, values, treat = NULL, cont_var_splits, split_type,
                       max_depth, min_samples, metric, int_weight_method,
                       erf_parent = NULL, trim = 1, d = 0, cluster = NULL) {
  # The ERF stored on each node is a thin <erf>: enough to plot the node's own
  # dose-response with plot.erf(), and nothing more. `erf_parent` itself keeps travelling
  # down the recursion as the full ERF function, because the split objective needs it.
  node_erf <- function() {
    thin_erf(erf_parent, values = values, treat = treat, split = split,
             bw = bw, size = sum(split), path = path)
  }

  if (depth >= max_depth) {
    # Maximum depth reached
    return(list(type = "leaf", path = path,
                size = sum(split), reason = "max_depth",
                erf = node_erf()))
  }

  if (sum(split) < min_samples * 2) {
    # Insufficient samples for splitting
    return(list(type = "leaf", path = path,
                size = sum(split), reason = "min_samples",
                erf = node_erf()))
  }

  # Make splitters from existing variables
  splitters <- do.call("c", lapply(names(M), function(v) {
    make_splitters(v, ss(M, split), cont_var_splits, split_type)
  }))

  if (is_null(splitters)) {
    # No effect modifiers available
    return(list(type = "leaf", path = path,
                size = sum(split), reason = "no_modifers",
                erf = node_erf()))
  }

  metrics <- alloc(NA_real_, length(splitters))

  for (i in seq_along(splitters)) {
    sp <- apply_splitter(splitters[[i]], split, M)

    if (sum(sp$splitL) < min_samples || sum(sp$splitR) < min_samples) {
      # Insufficient samples in subgroups
      next
    }

    # Splitting on modifier. A candidate that cannot be scored -- a degenerate subgroup
    # whose local-linear fit or integral weights fail -- costs only that candidate: its
    # metric stays NA and the search moves on, rather than aborting the whole tree. If
    # every candidate fails, the `allNA()` check below turns the node into a leaf.
    metrics[i] <- tryCatch({
      erfL <- estimate_erf_internal(fit, values = values,
                                    subset = sp$splitL,
                                    bw = bw, vcov = FALSE, verbose = FALSE, trim = trim,
                                    d = d)

      erfR <- estimate_erf_internal(fit, values = values,
                                    subset = sp$splitR,
                                    bw = bw, vcov = FALSE, verbose = FALSE, trim = trim,
                                    d = d)

      int_weights <- get_int_weights(erfL, erfR,
                                     values,
                                     int_weight_method)

      m <- get_metric(erfL, erfR,
                      values = values,
                      metric = metric,
                      int_weights = int_weights,
                      erf_parent = erf_parent)

      # A non-finite objective is not a usable candidate either.
      if (length(m) != 1L || !is.finite(m)) NA_real_ else m
    },
    error = function(e) NA_real_)
  }

  if (allNA(metrics)) {
    # No valid effect modifiers available
    return(list(type = "leaf", path = path,
                size = sum(split), reason = "no_valid_modifers",
                erf = node_erf()))
  }

  # Get best split
  best_i <- which(!is.na(metrics) & metrics == max(metrics, na.rm = TRUE))[1L]

  best_split <- apply_splitter(splitters[[best_i]], split, M)

  best_erfL <- estimate_erf_internal(fit, values = values,
                                     subset = best_split$splitL,
                                     bw = bw, vcov = TRUE, verbose = FALSE, trim = trim,
                                     d = d,
                                     cluster = subset_cluster(cluster, best_split$splitL))

  best_erfR <- estimate_erf_internal(fit, values = values,
                                     subset = best_split$splitR,
                                     bw = bw, vcov = TRUE, verbose = FALSE, trim = trim,
                                     d = d,
                                     cluster = subset_cluster(cluster, best_split$splitR))

  # Build tree recursively
  treeL <- Recall(M = M,
                  split = best_split$splitL,
                  depth = depth + 1L,
                  path = c(path, splitter_labels(splitters[[best_i]])[1L]),
                  fit = fit, bw = bw, values = values, treat = treat,
                  max_depth = max_depth, min_samples = min_samples,
                  metric = metric,
                  int_weight_method = int_weight_method,
                  cont_var_splits = cont_var_splits,
                  split_type = split_type,
                  erf_parent = best_erfL, trim = trim, d = d, cluster = cluster)

  treeR <- Recall(M = M,
                  split = best_split$splitR,
                  depth = depth + 1L,
                  path = c(path, splitter_labels(splitters[[best_i]])[2L]),
                  fit = fit, bw = bw, values = values, treat = treat,
                  max_depth = max_depth, min_samples = min_samples,
                  metric = metric,
                  int_weight_method = int_weight_method,
                  cont_var_splits = cont_var_splits,
                  split_type = split_type,
                  erf_parent = best_erfR, trim = trim, d = d, cluster = cluster)

  # Create internal node
  list(
    type = "internal",
    path = path,
    size = sum(split),
    split_var = splitter_variable(splitters[[best_i]]),
    split_levels = splitter_labels(splitters[[best_i]]),
    split_breaks = splitter_breaks(splitters[[best_i]]),
    split_levels_in = splitter_levels_in(splitters[[best_i]]),
    left = treeL,
    right = treeR,
    metric = metrics[best_i],
    erf = node_erf()
  )
}

# A "thin" <erf> for a single node of a tree: everything plot.erf() needs and nothing
# else. The node's ERF is re-expressed as its estimates and covariance at the anchor
# values, so the nuisance matrices that the full ERF's function carries in its
# environment are not retained once per node. Evaluated at the anchor values it returns
# exactly what the full ERF would.
thin_erf <- function(erf_fun, values, treat, split, bw = NULL, size = NULL, path = NULL) {
  if (!is.function(erf_fun)) {
    return(NULL)
  }

  est <- erf_fun(values, vcov = TRUE)

  node_treat <- treat[split]
  attr(node_treat, "varname") <- attr(treat, "varname")

  out <- list(erf = erf_interpolator(values, est[["est"]], est[["vcov"]]),
              values = values,
              treat = node_treat,
              bw = bw,
              size = size,
              path = path)

  class(out) <- c("erf_node", "erf")

  out
}

# The ERF function of a node, whether the node stores a thin <erf> or (as the
# development honest-tree builders do) the ERF function itself.
node_erf_fun <- function(node) {
  erf <- node[["erf"]]

  if (is.function(erf)) erf else erf[["erf"]]
}

apply_splitter <- function(splitter, current_split, data) {
  new_split <- splitter(data)

  list(splitL = current_split & new_split,
       splitR = current_split & !new_split)
}

make_cat_splitter <- function(levels, index_in, v) {
  splitter <- function(data) {
    data[[v]] %in% levels[index_in]
  }

  attr(splitter, "variable") <- v
  attr(splitter, "labels") <- make_labels(levels[index_in], levels[-index_in], v)

  # The level values sent to the left kid, recorded so the split is available
  # numerically and not only through its label.
  attr(splitter, "levels_in") <- levels[index_in]

  splitter
}

make_cont_splitter <- function(levels, index, unique, v) {
  splitter <- function(data) {
    # A unit missing this modifier goes to the right branch, which is what `%in%`
    # already does for a categorical modifier. Returning `NA` instead would carry
    # through `apply_splitter()` into the node sizes and abort the search.
    m <- data[[v]]
    !is.na(m) & m < levels[index]
  }

  attr(splitter, "variable") <- v
  attr(splitter, "labels") <- c(sprintf("%s %s %s",
                                        v,
                                        if (unique && index == 2L) "=" else "<",
                                        if (unique && index == 2L) levels[1L] else levels[index]),
                                sprintf("%s %s %s",
                                        v,
                                        if (unique && index == length(levels)) "=" else "\u2265",
                                        levels[index]))

  # The numeric cut point, recorded so the split is available numerically and not only
  # through its label. The comparison is `<`.
  attr(splitter, "breaks") <- levels[index]

  splitter
}

make_splitters <- function(v, data, num_splits, split_type) {
  m <- data[[v]]

  if (is.logical(m)) {
    m <- ordered(m)
  }

  if (is.ordered(m)) {
    m <- droplevels(m)

    u <- levels(m)

    if (length(u) == 1L) {
      return(NULL)
    }

    out <- lapply(seq_along(u)[-1L], function(i) {
      make_cat_splitter(u, seq_len(i - 1L), v)
    })

    return(out)
  }

  if (is.character(m)) {
    m <- factor(m)
  }

  if (is.factor(m)) {
    m <- droplevels(m)

    u <- levels(m)

    if (length(u) == 1L) {
      return(NULL)
    }

    combo_grid <- do.call(function(...) {expand.grid(..., KEEP.OUT.ATTRS = FALSE)},
                          setNames(rep.int(list(c(TRUE, FALSE)), length(u)), u)) |>
      as.matrix()

    combo_grid <- combo_grid[!(rowSums(combo_grid) %in% c(0, length(u))), , drop = FALSE]

    out <- lapply(seq_len(nrow(combo_grid)), function(i) {
      make_cat_splitter(u, which(combo_grid[i,]), v)
    })

    return(out)
  }

  u <- sort(unique(m))

  if (length(u) == 1L) {
    return(NULL)
  }

  if (length(u) <= num_splits) {
    out <- lapply(seq_along(u)[-1L], function(i) {
      make_cont_splitter(u, i, TRUE, v)
    })

    return(out)
  }

  if (split_type == "quantile") {
    # A modifier with a concentrated mass (a percentage that is zero for most units,
    # say) can put the same value at several quantiles, and `cut()` rejects duplicated
    # breaks. Dropping the duplicates leaves fewer bins than `num_splits`, which is the
    # most the data supports.
    breaks <- unique(c(-Inf,
                       quantile(m, seq_len(num_splits - 1L) / num_splits, na.rm = TRUE),
                       Inf))
  }
  else {
    breaks <- num_splits
  }

  b <- cut(m, breaks = breaks, right = FALSE) |>
    droplevels() |>
    unclass()

  # The cut points are midpoints between neighboring bins, so they are read off the
  # observed values only; a unit missing this modifier contributes to no bin.
  u <- vapply(seq_len(max(b, na.rm = TRUE) - 1L), function(i) {
    round_midpoint(max(m[which(b == i)]), min(m[which(b == i + 1L)]))
  }, numeric(1L))

  lapply(seq_along(u), function(i) {
    make_cont_splitter(u, i, FALSE, v)
  })
}

splitter_variable <- function(splitter) {
  attr(splitter, "variable")
}

splitter_labels <- function(splitter) {
  attr(splitter, "labels")
}

splitter_breaks <- function(splitter) {
  attr(splitter, "breaks")
}

splitter_levels_in <- function(splitter) {
  attr(splitter, "levels_in")
}

get_metric <- function(erfL, erfR, values, metric, int_weights, erf_parent = NULL) {
  estL <- erfL(values, vcov = FALSE)$est
  estR <- erfR(values, vcov = FALSE)$est

  integrate_erf <- function(x_vals, y_vals, weights) {
    weights_adjusted <- (weights[-1L] + weights[-length(weights)]) / 2

    sum(diff1(x_vals) * (y_vals[-1L] + y_vals[-length(y_vals)]) * weights_adjusted) / 2
  }

  switch(metric,
         "AUC" = {
           AUC_L <- integrate_erf(values, estL, int_weights)
           AUC_R <- integrate_erf(values, estR, int_weights)

           abs(AUC_L - AUC_R)
         },
         "L1" = {
           integrate_erf(values, abs(estL - estR), int_weights)
         },
         "L2" = {
           integrate_erf(values, (estL - estR)^2, int_weights)
         },
         "fANOVA" = {
           est_parent <- erf_parent(values, vcov = FALSE)$est

           L2_L <- integrate_erf(values, (estL - est_parent)^2, int_weights)
           L2_R <- integrate_erf(values, (estR - est_parent)^2, int_weights)

           # Each child's deviation from the parent is weighted by that child's share of
           # the node, so a small subgroup that happens to deviate wildly does not score
           # like a large one. Weighted counts rather than row counts, to match the
           # weighted estimator; these equal the row counts when unweighted. The split is
           # exhaustive, so the two shares sum to the parent.
           nL <- sum(attr(erfL, "res")[["weights"]])
           nR <- sum(attr(erfR, "res")[["weights"]])

           (nL * L2_L + nR * L2_R) / (nL + nR)
         })
}

get_int_weights <- function(erfL, erfR, values, int_weight_method) {
  if (int_weight_method == "none") {
    return(alloc(1.0, length(values)))
  }

  AL <- attr(erfL, "res")[["A"]]
  AR <- attr(erfR, "res")[["A"]]

  weightsL <- attr(erfL, "res")[["weights"]]
  weightsR <- attr(erfR, "res")[["weights"]]

  # Interpolate a density onto the anchor values. Anything outside the range the density
  # is defined over gets 0, not a flat extrapolation: the exposure density really is ~0
  # where the subgroup has no support, and w() below is a common-support weight, so 0 is
  # what makes it one.
  interp_dens <- function(x, y) {
    out <- approxfun(x = x, y = y, rule = 1)(values)
    out[is.na(out) | out < 0] <- 0
    out
  }

  if (int_weight_method == "density") {
    # Calculate densities
    densL <- density(AL,
                     from = min(values),
                     to = max(values),
                     weights = weightsL / sum(weightsL))
    densR <- density(AR,
                     from = min(values),
                     to = max(values),
                     weights = weightsR / sum(weightsR))

    # Extract values
    piL <- interp_dens(densL$x, densL$y)
    piR <- interp_dens(densR$x, densR$y)
  }
  else if (int_weight_method == "gps") {
    # Average the conditional exposure density over the units in each subgroup. This is a
    # column mean of `pi_hat_mat`, so it is indexed by the nuisance exposure grid
    # `A_vals` -- not by the units' own exposures -- and must be interpolated as such.
    A_vals <- attr(erfL, "res")[["A_vals"]]

    phatL <- fmean(attr(erfL, "res")[["pi_hat_mat"]], w = weightsL)
    phatR <- fmean(attr(erfR, "res")[["pi_hat_mat"]], w = weightsR)

    # Interpolate to fill in the grid
    piL <- interp_dens(A_vals, phatL)
    piR <- interp_dens(A_vals, phatR)
  }

  incr <- mean(diff1(values))

  # Normalize densities
  piL <- piL / sum(piL * incr)
  piR <- piR / sum(piR * incr)

  # Calculate weight function
  w_unnormalized <- piL * piR
  normalizing_constant <- sum(w_unnormalized * incr)

  # Two subgroups can have effectively disjoint exposure support, leaving no anchor value
  # where both have mass and a normalizing constant of 0. Weight the anchors uniformly
  # rather than returning NaN, which would silently drop the candidate split.
  if (!is.finite(normalizing_constant) || normalizing_constant <= 0) {
    return(alloc(1.0, length(values)))
  }

  w_unnormalized / normalizing_constant
}

# Helper function to apply a tree's learned rules on a new ERF.
apply_tree <- function(tree_grow, erf_est) {
  check_full_erf(erf_est)
  
  modifiers <- stats::reformulate(names(tree_grow$M))
  M <- process_modifiers(modifiers, erf_est)
  
  values <- erf_est$values
  d <- if (is.null(erf_est$d)) 0L else erf_est$d
  trim <- if (is.null(erf_est$trim)) 1 else erf_est$trim
  
  # Reconstruct the exact stored rule.
  get_splitter <- function(node) {
    v <- node$split_var

    if (!is.null(node$split_levels_in)) {
      levels_in <- node$split_levels_in
      
      function(data) {
        data[[v]] %in% levels_in
      }
    } else if (!is.null(node$split_breaks)) {
      cutoff <- node$split_breaks
     
      function(data) {
        x <- data[[v]]
        !is.na(x) & x < cutoff
      }
    } else {
      stop("Node has no stored split rule.")
    }
  }
  
  visit <- function(node, split) {
    
    ef <- estimate_erf_internal(
      fit = erf_est$fit,
      values = values,
      subset = split,
      bw = erf_est$bw,
      vcov = TRUE,
      verbose = FALSE,
      trim = trim,
      d = d
    )
    
    node$size_grow <- node$size
    node$size <- sum(split)
    
    node$erf <- thin_erf(
      erf_fun = ef,
      values = values,
      treat = erf_est$treat,
      split = split,
      bw = attr(ef, "bw"),
      size = node$size,
      path = node$path
    )
    
    if (node$type == "internal") {
      sp <- apply_splitter(
        splitter = get_splitter(node),
        current_split = split,
        data = M
      )
      
      node$left  <- visit(node$left, sp$splitL)
      node$right <- visit(node$right, sp$splitR)
    }
    
    node
  }
  
  out <- tree_grow
  out$tree <- visit(tree_grow$tree, rep(TRUE, nrow(M)))
  
  # Replace discovery-sample metadata with estimation-sample metadata.
  out$M <- M
  out$method <- erf_est$method
  out$values <- values
  out$treat <- erf_est$treat
  out$treat_name <- attr(erf_est$treat, "varname")
  out$pseudo_Y <- attr(erf_est$erf, "res")[["pseudo_Y"]]
  out$weights <- erf_est$weights
  
  out
}

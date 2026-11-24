#' Grow an Honest ERF Tree
#'
#' @description
#' `honest_erf_tree()` discovers subgroups with differing ERFs on one part of the sample
#' and estimates the subgroup ERFs on the other, so that the confidence bands at the
#' nodes of the tree are not narrowed by the search that produced it. This involves
#' splitting the sample into a discovery sample and an inference sample with
#' [split_sample()], estimating an ERF on each with [estimate_erf()], growing the tree on
#' the discovery ERF with [rpart_erf()], and re-estimating the ERF at every node of that
#' tree on the inference sample with [refit_erf_tree()].
#'
#' @inheritParams estimate_erf
#' @inheritParams rpart_erf
#' @inheritParams split_sample
#' @param modifiers a one-sided formula containing the potential effect modifiers, e.g.,
#'   `~ M1 + M2`. The modifiers must be among the covariates in `formula`, since they are
#'   read from the `data` component of the ERF estimated on each sample. If `NULL`, the
#'   default, all covariates are searched. A unit missing the modifier on which a split is made
#'   is sent to the right branch of that split.
#' @param cluster an optional one-sided formula (e.g., `~county`), vector, or one-column
#'   data frame giving the clustering variable, with one value per row of `data`. When
#'   supplied, whole clusters are assigned to the discovery or the inference sample
#'   together, the cross-fitting folds within each sample are formed over whole clusters,
#'   and the variance of every node ERF is cluster-robust. Default is `NULL` for no
#'   clustering.
#' @param ... additional arguments passed to [estimate_erf()].
#'
#' @return
#' An `<honest_erf_tree>` object, which is an `<erf_tree>` object (see [rpart_erf()])
#' whose node ERFs are estimated on the inference sample, with the additional components
#' described at [refit_erf_tree()] and two more:
#' \item{`split`}{a factor with one entry per row of `data` and levels `"discovery"` and
#' `"inference"`, giving the sample each row was assigned to}
#' \item{`call`}{the original function call}
#'
#' @details
#' A tree grown and evaluated on the same data reports bands at its nodes that do not
#' account for the search over candidate splits that selected those nodes, because the
#' split that maximizes the difference between two subgroup ERFs also tends to be the one
#' whose estimated difference is largest by chance. Honest sample splitting (Athey &
#' Imbens, 2016) separates the two tasks. The discovery sample determines the partition,
#' and the inference sample, which played no part in choosing it, supplies the subgroup
#' curves and their pointwise confidence bands. Conditional on the tree, those bands then
#' have the same justification as the bands of an ERF estimated by [estimate_erf()] on a
#' subgroup chosen in advance.
#'
#' The procedure follows the two phases of the algorithm in the methods paper the package
#' implements. The sample is split at random, by cluster when `cluster` is supplied. On
#' the discovery sample the nuisance functionals are estimated, cross-fit when `crossfit`
#' is greater than 1, and the tree is grown by recursive partitioning, at each node
#' recomputing the subgroup pseudo-outcomes of every candidate split and maximizing the
#' objective named by `metric`. On the inference sample the nuisance functionals are
#' estimated afresh, and the pseudo-outcomes and the ERF are recomputed within each node
#' of the tree that was grown. Nothing estimated on one sample is carried into the other,
#' so the two ERFs do not share a bandwidth unless `bw` is supplied as a number: when it
#' is `NULL`, each ERF selects its own by leave-one-out cross-validation on its own sample.
#'
#' The price of honesty is sample size. Each ERF is estimated on a fraction of the data,
#' so the discovery sample has less power to find a split than the full sample would, and
#' the bands at the nodes are wider than those of a tree grown on everything. `frac`
#' controls the division; a larger discovery sample finds splits more reliably, and a
#' larger inference sample estimates the node curves more precisely. The exposure grid at
#' which the node ERFs are held is the `values` component of the inference ERF, which is
#' determined by `range` and `n` on that sample.
#'
#' The bands are pointwise, cover the curve of each node on its own rather than the whole
#' tree simultaneously, and are valid under the conditions of [estimate_erf()] when the
#' observations are independent. When they are clustered, assigning whole clusters to one
#' sample keeps the two samples independent of one another and the node variances
#' cluster-robust, but dependence between clusters, such as spatial dependence between
#' neighboring counties, is not accounted for.
#'
#' To carry out the two phases separately, as when the inference ERF should use a
#' different nuisance method or the discovery tree is to be inspected before any curve is
#' estimated on the inference sample, call [split_sample()], [estimate_erf()] on each
#' sample, [rpart_erf()] on the discovery ERF, and [refit_erf_tree()] with the inference
#' ERF. See Examples.
#'
#' @seealso [refit_erf_tree()] for the inference phase on its own; [split_sample()];
#' [rpart_erf()] for growing a tree without sample splitting; [summary.erf_tree()],
#' [plot.erf_tree()], and [`[[.erf_tree`] for reading the resulting tree.
#'
#' @references
#' Athey, S., & Imbens, G. (2016). Recursive partitioning for heterogeneous causal
#' effects. *Proceedings of the National Academy of Sciences*, 113(27), 7353–7360.
#' \doi{10.1073/pnas.1510489113}
#'
#' Kennedy, E. H., Ma, Z., McHugh, M. D., & Small, D. S. (2017). Non-parametric methods for doubly robust estimation of continuous treatment effects. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 79(4), 1229–1245. \doi{10.1111/rssb.12212}
#'
#' @examples
#' # Simulated data with an effect modifier: the exposure-
#' # response slope differs by the binary modifier M.
#' set.seed(1234)
#' n <- 1000
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' M <- factor(rbinom(n, 1, 0.5))
#' A <- 1 + 0.5 * X1 + rnorm(n)
#' Y <- A * (1 + 1.5 * (M == "1")) + X1 + X2 + rnorm(n)
#' d <- data.frame(Y, A, X1, X2, M)
#'
#' # Half the sample grows the tree; the other half estimates its ERFs
#' ht <- honest_erf_tree(Y ~ A | X1 + X2 + M, data = d,
#'                       modifiers = ~ M + X1,
#'                       max_depth = 2, min_samples = 50,
#'                       n = 25L)
#' ht
#' summary(ht)
#'
#' # The ERF at a node, estimated on the inference sample
#' plot(ht[[2]])
#'
#' \donttest{
#' # The same two phases, step by step
#' set.seed(1234)
#' s <- split_sample(d)
#'
#' erf_dis <- estimate_erf(Y ~ A | X1 + X2 + M,
#'                         data = d[s == "discovery", ],
#'                         n = 25L)
#'
#' tr <- rpart_erf(erf_dis, ~ M + X1, max_depth = 2,
#'                 min_samples = 50)
#'
#' erf_inf <- estimate_erf(Y ~ A | X1 + X2 + M,
#'                         data = d[s == "inference", ],
#'                         n = 25L)
#'
#' ht2 <- refit_erf_tree(tr, erf_inf)
#'
#' # The tree the search found, and the whole tree with honest node ERFs
#' summary(ht2$discovery)
#' plot(ht2)
#' }
#'
#' @export
honest_erf_tree <- function(formula, data, modifiers = NULL, method = method_glm(),
                            weights = NULL, bw = NULL, crossfit = FALSE, cluster = NULL,
                            n = 51L, range = 1, n_grid = 100L, trim = 0.99, d = 0,
                            metric = "L2", min_samples = 30, max_depth = 3,
                            int_weight_method = "density", cont_var_splits = 5L,
                            split_type = "quantile", frac = 0.5, verbose = FALSE, ...) {
  call <- match.call()

  arg::arg_data.frame(data)

  # The tree settings are checked before anything is fit, so a bad one does not surface
  # only after two nuisance fits.
  metric <- process_metric(metric)
  int_weight_method <- process_int_weight_method(int_weight_method)
  arg::arg_count(min_samples)
  arg::arg_gt(min_samples, 1)
  arg::arg_count(max_depth)
  arg::arg_count(cont_var_splits)
  split_type <- arg::match_arg(split_type, c("quantile", "range"))

  weights <- process_weights(weights, data)
  cluster <- resolve_cluster(cluster, data)

  split <- split_sample(data, cluster = cluster, frac = frac)

  dis <- split == "discovery"

  if (verbose) {
    message(sprintf("Discovery sample: %d rows; inference sample: %d rows",
                    sum(dis), sum(!dis)))
  }

  fit_half <- function(rows) {
    estimate_erf(formula, data = data[rows, , drop = FALSE], method = method,
                 weights = weights[rows], bw = bw, crossfit = crossfit,
                 cluster = if (is_null(cluster)) NULL else cluster[rows, , drop = FALSE],
                 n = n, range = range, n_grid = n_grid, trim = trim,
                 verbose = verbose, d = d, ...)
  }

  if (verbose) {
    message("Discovery phase: estimating the ERF and growing the tree...")
  }

  erf_dis <- fit_half(dis)

  tree <- rpart_erf(erf_dis, modifiers = modifiers, metric = metric,
                    min_samples = min_samples, max_depth = max_depth,
                    int_weight_method = int_weight_method,
                    cont_var_splits = cont_var_splits, split_type = split_type)

  if (verbose) {
    message("Inference phase: estimating the ERF and the node ERFs...")
  }

  erf_inf <- fit_half(!dis)

  out <- refit_erf_tree(tree, erf_inf)

  out[["split"]] <- split
  out[["call"]] <- call

  out
}

#' Re-estimate the ERFs of a Tree on an Independent Sample
#'
#' @description
#' `refit_erf_tree()` applies the splits of an `<erf_tree>` to the sample an `<erf>` was
#' estimated on and estimates the ERF at every node of the tree there. When that sample
#' took no part in growing the tree, the result is an honest tree, whose node ERFs and
#' confidence bands are estimated on data the search never saw. This is the inference
#' phase of [honest_erf_tree()], available on its own.
#'
#' @param tree an `<erf_tree>` object; the output of a call to [rpart_erf()], grown on
#'   the discovery sample.
#' @param erf an `<erf>` object; the output of a call to [estimate_erf()] on the
#'   inference sample. Its `data` component must contain every modifier the tree splits
#'   on, and it must have been estimated at the same order `d` as the ERF the tree was
#'   grown from.
#'
#' @return
#' An `<honest_erf_tree>` object, which is an `<erf_tree>` object (see [rpart_erf()])
#' with the same tree structure as `tree` and the following differences. The `size` of
#' every node, the `erf` held at it, the `M`, `values`, `treat`, `pseudo_Y`, and
#' `weights` components, and the `metric` of every internal node all come from the
#' inference sample, and the following components are added:
#' \item{`discovery`}{the `<erf_tree>` supplied as `tree`, with its node ERFs estimated
#' on the discovery sample}
#' Each node additionally carries `size_discovery`, the number of discovery-sample units
#' at that node, and each internal node `metric_discovery`, the value of the objective
#' for its split on the discovery sample.
#'
#' @details
#' The tree is descended from the root. At each internal node, the units of the
#' inference sample are sent left or right by the recorded split: for a continuous
#' modifier, units below the cut point go left, and for a categorical one, units in the
#' recorded set of levels go left; a unit missing the modifier goes right, as when the
#' tree was grown. Within each node the pseudo-outcomes are recomputed from the nuisance
#' estimates stored in `erf`, so that the stabilizing density and the g-computation
#' component are those of the subgroup, and the ERF is fit by local linear regression at
#' the bandwidth of `erf`, with its trimming of the density ratio and its clustering. The
#' ERF at the root is the ERF of `erf` itself. The anchor values are those of `erf`.
#'
#' The objective that selected each split is evaluated again on the inference sample
#' from the two child ERFs, with the metric and the integral weighting the tree was grown
#' with, and stored as that node's `metric`; the value on the discovery sample, which is
#' what the search maximized, is kept as `metric_discovery`. The difference between the
#' two is a measure of how much the search inflated the apparent heterogeneity.
#'
#' A node that holds fewer than two units of the inference sample has no ERF, and an
#' error is thrown naming it; a node that holds fewer units than the `min_samples` the
#' tree was grown with produces a warning, since its band rests on fewer units than the
#' tree required of a leaf. Both are more likely with a deep tree or a small inference
#' sample.
#'
#' Supplying the very ERF the tree was grown from returns the tree with the node ERFs
#' [rpart_erf()] already computed, which makes the function a check on the recorded
#' splits but adds nothing to the inference.
#'
#' @seealso [honest_erf_tree()], which runs both phases; [split_sample()];
#' [rpart_erf()]; [summary.erf_tree()], [plot.erf_tree()], and [`[[.erf_tree`] for
#' reading the resulting tree.
#'
#' @references
#' Athey, S., & Imbens, G. (2016). Recursive partitioning for heterogeneous causal
#' effects. *Proceedings of the National Academy of Sciences*, 113(27), 7353–7360.
#' \doi{10.1073/pnas.1510489113}
#'
#' @examples
#' set.seed(1234)
#' n <- 1000
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' M <- factor(rbinom(n, 1, 0.5))
#' A <- 1 + 0.5 * X1 + rnorm(n)
#' Y <- A * (1 + 1.5 * (M == "1")) + X1 + X2 + rnorm(n)
#' d <- data.frame(Y, A, X1, X2, M)
#'
#' # Split the sample, and estimate an ERF on each half
#' s <- split_sample(d)
#'
#' erf_dis <- estimate_erf(Y ~ A | X1 + X2 + M, data = d[s == "discovery", ], n = 25L)
#' erf_inf <- estimate_erf(Y ~ A | X1 + X2 + M, data = d[s == "inference", ], n = 25L)
#'
#' # Grow the tree on the discovery half
#' tr <- rpart_erf(erf_dis, ~ M + X1, max_depth = 2, min_samples = 50)
#'
#' # Estimate its node ERFs on the inference half
#' ht <- refit_erf_tree(tr, erf_inf)
#' ht
#'
#' # Node sizes and the objective at the root, on each sample
#' summary(ht)
#' c(discovery = ht$tree$metric_discovery, inference = ht$tree$metric)
#'
#' @export
refit_erf_tree <- function(tree, erf) {
  arg::arg_is(tree, "erf_tree")
  arg::arg_is(erf, "erf")
  check_full_erf(erf)

  d_tree <- tree[["d"]] %or% 0L
  d <- erf[["d"]] %or% 0L

  if (d_tree != d) {
    arg::err("{.arg tree} was grown from an ERF estimated with {.code d = {d_tree}}, but {.arg erf} was estimated with {.code d = {d}}; both must be the ERF itself or both its derivative")
  }

  # The modifiers are resolved against the inference ERF the same way rpart_erf() resolved
  # them against the discovery ERF, so they must be columns of `erf$data`.
  M <- process_modifiers(tree[["modifiers"]], erf)

  split_vars <- unique(tree_split_vars(tree[["tree"]]))

  missing_vars <- setdiff(split_vars, names(M))

  if (!is_null(missing_vars)) {
    arg::err("the tree splits on {.var {missing_vars}}, which {?is/are} not among the modifiers available in {.arg erf}")
  }

  if (nrow(M) != length(erf[["treat"]])) {
    arg::err("the modifiers in {.arg erf} do not have one row per unit the ERF was estimated on")
  }

  fit <- erf[["fit"]]
  values <- erf[["values"]]
  treat <- erf[["treat"]]
  bw <- erf[["bw"]]
  trim <- erf[["trim"]] %or% 1
  cluster <- erf[["cluster"]]

  metric <- tree[["info"]][["metric"]]
  int_weight_method <- tree[["info"]][["int_weight_method"]]
  min_samples <- tree[["info"]][["min_samples"]] %or% 0

  small <- character(0)

  node_erf_fun <- function(split, path) {
    n_node <- sum(split)
    label <- path_label(path) %or% "root"

    if (n_node < 2L) {
      arg::err(c("Node {.val {label}} holds {n_node} unit{?s} of the inference sample, too few to estimate its ERF.",
                 i = "A deeper tree or a smaller inference sample leaves fewer units at each node; use a larger {.arg frac}, a shallower tree, or a larger {.arg min_samples}."))
    }

    if (n_node < min_samples) {
      small <<- c(small, label)
    }

    tryCatch({
      estimate_erf_internal(fit, values = values, subset = split, bw = bw, vcov = TRUE,
                            cluster = subset_cluster(cluster, split), verbose = FALSE,
                            trim = trim, d = d)
    },
    error = function(e) {
      arg::err(c("The ERF could not be estimated at node {.val {label}} ({n_node} unit{?s} of the inference sample).",
                 x = conditionMessage(e)))
    })
  }

  # Descend the tree, sending the inference units down the recorded splits. The ERF
  # function of each node is passed in by its parent, which needs both children's ERFs to
  # evaluate the split objective on the inference sample before recursing.
  refit_node <- function(node, split, erf_fun) {
    size <- sum(split)

    thin <- thin_erf(erf_fun, values = values, treat = treat, split = split,
                     bw = bw, size = size, path = node[["path"]])

    if (identical(node[["type"]], "leaf")) {
      return(list(type = "leaf",
                  path = node[["path"]],
                  size = size,
                  size_discovery = node[["size"]],
                  reason = node[["reason"]],
                  erf = thin))
    }

    left <- split & node_goes_left(node, M)
    right <- split & !node_goes_left(node, M)

    erfL <- node_erf_fun(left, node[["left"]][["path"]])
    erfR <- node_erf_fun(right, node[["right"]][["path"]])

    # The objective on the inference sample. A candidate that cannot be scored here is
    # recorded as NA rather than aborting the refit, as in the search itself.
    met <- tryCatch({
      int_weights <- get_int_weights(erfL, erfR, values, int_weight_method)

      m <- get_metric(erfL, erfR, values = values, metric = metric,
                      int_weights = int_weights, erf_parent = erf_fun)

      if (length(m) != 1L || !is.finite(m)) NA_real_ else m
    },
    error = function(e) NA_real_)

    list(type = "internal",
         path = node[["path"]],
         size = size,
         size_discovery = node[["size"]],
         split_var = node[["split_var"]],
         split_levels = node[["split_levels"]],
         split_breaks = node[["split_breaks"]],
         split_levels_in = node[["split_levels_in"]],
         left = refit_node(node[["left"]], left, erfL),
         right = refit_node(node[["right"]], right, erfR),
         metric = met,
         metric_discovery = node[["metric"]],
         erf = thin)
  }

  # The root holds every inference unit, so its ERF is the inference ERF itself.
  root <- refit_node(tree[["tree"]], rep.int(TRUE, nrow(M)), erf[["erf"]])

  if (!is_null(small)) {
    arg::wrn(c("{length(small)} node{?s} of the tree hold{?s/} fewer than {min_samples} unit{?s} of the inference sample: {.val {small}}.",
               i = "The {cli::qty(length(small))} ERF{?s} there rest{?s/} on fewer units than the tree required of a leaf; the band{?s} may be unreliable."))
  }

  out <- list(tree = root,
              M = M,
              method = erf[["method"]],
              values = values,
              treat = treat,
              treat_name = attr(treat, "varname"),
              pseudo_Y = attr(erf[["erf"]], "res")[["pseudo_Y"]],
              weights = erf[["weights"]],
              modifiers = tree[["modifiers"]],
              d = d,
              info = tree[["info"]],
              discovery = tree)

  class(out) <- c("honest_erf_tree", "erf_tree")

  out
}

#' Split a Sample into Discovery and Inference Samples
#'
#' @description
#' `split_sample()` assigns each row of a dataset at random to a discovery sample or an
#' inference sample, keeping whole clusters together when a clustering variable is
#' supplied. The result is the sample split that [honest_erf_tree()] performs, for use
#' when its two phases are run separately.
#'
#' @param data a data frame; the dataset to split.
#' @param cluster an optional one-sided formula (e.g., `~county`), vector, or one-column
#'   data frame giving the clustering variable, with one value per row of `data`. When
#'   supplied, whole clusters are assigned to one sample or the other. Default is `NULL`
#'   to assign rows individually.
#' @param frac `numeric`; the fraction of the rows, or of the clusters when `cluster` is
#'   supplied, assigned to the discovery sample. Default is .5 to split the sample in
#'   half.
#'
#' @return
#' A factor with one entry per row of `data` and levels `"discovery"` and `"inference"`.
#'
#' @details
#' Without `cluster`, `floor(frac * nrow(data))` rows are sampled without replacement
#' into the discovery sample and the rest form the inference sample. With `cluster`, the
#' same is done to the distinct clusters, so that every row of a cluster lands in the
#' same sample and the fraction of rows in the discovery sample matches `frac` only
#' approximately. The split is meant to make the two samples independent, and repeated
#' observations of the same unit would otherwise let the inference sample share
#' information with the discovery sample that chose the tree.
#'
#' The assignment is random; call [set.seed()] first to make it reproducible.
#'
#' @seealso [honest_erf_tree()], which calls `split_sample()`; [refit_erf_tree()]
#'
#' @examples
#' data("nhanes3lead", package = "adrftools")
#'
#' set.seed(1234)
#' s <- split_sample(nhanes3lead)
#' table(s)
#'
#' # Clustered data: whole clusters stay together, so the split is at
#' # the cluster level
#' nhanes3lead$id <- factor(rep(seq_len(50), length.out = nrow(nhanes3lead)))
#'
#' s_cl <- split_sample(nhanes3lead, cluster = ~id)
#' table(s_cl)
#' table(nhanes3lead$id, s_cl)[1:4, ]
#'
#' @export
split_sample <- function(data, cluster = NULL, frac = .5) {
  arg::arg_data.frame(data)
  arg::arg_number(frac)
  arg::arg_between(frac, c(0, 1), inclusive = FALSE)

  n <- nrow(data)

  cluster <- resolve_cluster(cluster, data)

  if (is_null(cluster)) {
    k <- floor(frac * n)

    if (k < 1L || k >= n) {
      arg::err("{.arg frac} leaves no rows in one of the two samples")
    }

    discovery <- seq_len(n) %in% sample.int(n, k)
  }
  else {
    cl <- as.character(cluster[[1L]])
    lev <- unique(cl)
    k <- floor(frac * length(lev))

    if (k < 1L || k >= length(lev)) {
      arg::err("{.arg frac} leaves no clusters in one of the two samples")
    }

    discovery <- cl %in% sample(lev, k)
  }

  factor(ifelse(discovery, "discovery", "inference"),
         levels = c("discovery", "inference"))
}

# Whether each unit of a modifier frame goes to the left child of an internal node, read
# off the split the node records: below the cut point for a continuous modifier, in the
# recorded set of levels for a categorical one. A unit missing the modifier goes right in
# both cases, as in the splitters that grew the tree.
node_goes_left <- function(node, M) {
  m <- M[[node[["split_var"]]]]

  if (!is_null(node[["split_breaks"]])) {
    return(!is.na(m) & m < node[["split_breaks"]])
  }

  m %in% node[["split_levels_in"]]
}

# The variables a tree splits on, in depth-first order.
tree_split_vars <- function(node) {
  if (identical(node[["type"]], "leaf")) {
    return(character(0))
  }

  c(node[["split_var"]], Recall(node[["left"]]), Recall(node[["right"]]))
}

# The clustering of an <erf>, restricted to the units in a node and processed as
# estimate_erf() processes it, so that the number of clusters the small-sample adjustment
# counts is the number present in the node.
subset_cluster <- function(cluster, split) {
  if (is_null(cluster)) {
    return(NULL)
  }

  cl <- ss(cluster, split)

  process_cluster(cl, cl)
}

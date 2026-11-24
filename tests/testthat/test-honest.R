# The honest tree is specified by provenance: every node ERF and every node size must come
# from the inference sample, and the structure must be the discovery tree's. These helpers
# walk a tree independently of the package's own utilities so the two sides cannot drift
# together.

walk_tree <- function(node) {
  out <- list(node)

  if (!is.null(node$left)) {
    out <- c(out, walk_tree(node$left))
  }

  if (!is.null(node$right)) {
    out <- c(out, walk_tree(node$right))
  }

  out
}

split_vars <- function(node) {
  if (identical(node$type, "leaf")) {
    return(character(0))
  }

  c(node$split_var, split_vars(node$left), split_vars(node$right))
}

fit_halves <- function(d, n = 21L, ...) {
  s <- split_sample(d, ...)

  list(split = s,
       dis = estimate_erf(Y ~ A | X1 + X2 + mod, data = d[s == "discovery", ], n = n),
       inf = estimate_erf(Y ~ A | X1 + X2 + mod, data = d[s == "inference", ], n = n))
}

test_that("split_sample assigns rows, or whole clusters, at the requested fraction", {
  d <- sim_erf_data(n = 400, clsize = 10)

  set.seed(1)
  s <- split_sample(d, frac = .3)

  expect_s3_class(s, "factor")
  expect_identical(levels(s), c("discovery", "inference"))
  expect_length(s, nrow(d))
  expect_identical(sum(s == "discovery"), 120L)

  # Clustered: no cluster straddles the two samples, and the fraction applies to the
  # clusters (40 of them) rather than to the rows.
  set.seed(1)
  s_cl <- split_sample(d, cluster = ~cl, frac = .5)

  spread <- tapply(s_cl, d$cl, function(z) length(unique(z)))
  expect_true(all(spread == 1L))
  expect_identical(length(unique(d$cl[s_cl == "discovery"])), 20L)

  # Reproducible under a seed, and the same for the three ways of supplying the clusters.
  set.seed(1)
  expect_identical(split_sample(d, cluster = d$cl), s_cl)
  set.seed(1)
  expect_identical(split_sample(d, cluster = d["cl"]), s_cl)

  # A fraction that empties one sample is refused rather than returned.
  expect_error(split_sample(d, frac = 0))
  expect_error(split_sample(d, frac = 1))
  expect_error(split_sample(d, frac = 1 / 1000), "no rows")
  expect_error(split_sample(d, cluster = ~cl, frac = 1 / 100), "no clusters")
  expect_error(split_sample(as.matrix(d[1:2])))
})

test_that("refitting a tree on the ERF it was grown from reproduces rpart_erf", {
  d <- sim_erf_data(n = 500)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 50)

  rf <- refit_erf_tree(tr, e)

  expect_s3_class(rf, "honest_erf_tree")
  expect_s3_class(rf, "erf_tree")

  nodes_tr <- walk_tree(tr$tree)
  nodes_rf <- walk_tree(rf$tree)

  expect_length(nodes_rf, length(nodes_tr))

  for (i in seq_along(nodes_tr)) {
    expect_identical(nodes_rf[[i]]$type, nodes_tr[[i]]$type)
    expect_identical(nodes_rf[[i]]$path, nodes_tr[[i]]$path)
    expect_identical(nodes_rf[[i]]$size, nodes_tr[[i]]$size)
    expect_identical(nodes_rf[[i]]$size_discovery, nodes_tr[[i]]$size)

    # The same subset, the same nuisance, the same bandwidth: the same curve.
    est_tr <- nodes_tr[[i]]$erf$erf(tr$values, vcov = TRUE)
    est_rf <- nodes_rf[[i]]$erf$erf(rf$values, vcov = TRUE)
    expect_equal(est_rf$est, est_tr$est)
    expect_equal(est_rf$vcov, est_tr$vcov)
  }

  expect_identical(split_vars(rf$tree), split_vars(tr$tree))

  # The objective on the "inference" sample is the one the search maximized.
  expect_equal(rf$tree$metric, rf$tree$metric_discovery)
  expect_equal(rf$tree$metric, tr$tree$metric)
})

test_that("node sizes, membership, and ERFs come from the inference sample alone", {
  d <- sim_erf_data(n = 800)

  set.seed(2)
  h <- fit_halves(d)

  tr <- rpart_erf(h$dis, ~ mod + X1, max_depth = 2, min_samples = 50)
  ht <- refit_erf_tree(tr, h$inf)

  n_inf <- sum(h$split == "inference")

  # The root holds every inference unit and the children partition it.
  expect_identical(ht$tree$size, n_inf)
  expect_identical(ht$tree$size, ht$tree$left$size + ht$tree$right$size)
  expect_identical(nrow(ht$M), n_inf)
  expect_identical(length(ht$treat), n_inf)

  # The structure is the discovery tree's, and the discovery sizes are kept beside the
  # inference sizes.
  expect_identical(split_vars(ht$tree), split_vars(tr$tree))
  expect_identical(ht$tree$size_discovery, tr$tree$size)
  expect_identical(ht$tree$size_discovery, sum(h$split == "discovery"))
  expect_identical(ht$discovery$tree$size, tr$tree$size)

  # Membership at the root split, applied by hand to the inference data from the recorded
  # cut, matches what the left node holds.
  M_inf <- h$inf$data
  root <- ht$tree
  left <- {
    if (!is.null(root$split_breaks)) !is.na(M_inf[[root$split_var]]) & M_inf[[root$split_var]] < root$split_breaks
    else M_inf[[root$split_var]] %in% root$split_levels_in
  }
  expect_identical(root$left$size, sum(left))
  expect_equal(sort(as.numeric(root$left$erf$treat)), sort(as.numeric(h$inf$treat[left])))

  # The root ERF is the inference ERF, not the discovery one. The two ERFs have their
  # own anchor ranges, so the discovery one is evaluated only where both are defined.
  expect_equal(ht[[1L]]$erf(ht$values)$est, h$inf$erf(h$inf$values)$est)

  inside <- seq(max(min(h$dis$values), min(ht$values)),
                min(max(h$dis$values), max(ht$values)), length.out = 7L)
  expect_false(isTRUE(all.equal(ht[[1L]]$erf(inside, vcov = FALSE)$est,
                                h$dis$erf(inside, vcov = FALSE)$est)))

  # The anchor values are the inference ERF's.
  expect_identical(ht$values, h$inf$values)
})

test_that("honest node ERFs recover the subgroup slopes", {
  # mod = 0 has slope 1 and mod = 1 has slope 2.5; the tree should split on mod and the
  # inference-sample curves at its two children should carry those slopes.
  d <- sim_erf_data(n = 2000, seed = 11)

  set.seed(3)
  ht <- honest_erf_tree(Y ~ A | X1 + X2 + mod, data = d, modifiers = ~mod,
                        max_depth = 1, min_samples = 50, n = 21L)

  expect_identical(ht$tree$split_var, "mod")

  slope <- function(node) {
    est <- node$erf$erf(c(0, 2), vcov = FALSE)$est
    (est[2L] - est[1L]) / 2
  }

  left_is_0 <- "0" %in% ht$tree$split_levels_in
  s0 <- slope(if (left_is_0) ht$tree$left else ht$tree$right)
  s1 <- slope(if (left_is_0) ht$tree$right else ht$tree$left)

  expect_equal(s0, 1, tolerance = .25)
  expect_equal(s1, 2.5, tolerance = .15)
})

test_that("the split objective is evaluated again on the inference sample", {
  d <- sim_erf_data(n = 800)

  set.seed(4)
  h <- fit_halves(d)

  tr <- rpart_erf(h$dis, ~ mod + X1, max_depth = 2, min_samples = 50)
  ht <- refit_erf_tree(tr, h$inf)

  internal <- Filter(function(nd) identical(nd$type, "internal"), walk_tree(ht$tree))

  for (nd in internal) {
    expect_type(nd$metric, "double")
    expect_true(is.finite(nd$metric))
    expect_true(is.finite(nd$metric_discovery))
    expect_false(isTRUE(all.equal(nd$metric, nd$metric_discovery)))
  }

  # The discovery value is the one the search recorded.
  expect_equal(ht$tree$metric_discovery, tr$tree$metric)
})

test_that("honest_erf_tree splits by cluster and gives cluster-robust node bands", {
  d <- sim_erf_data(n = 600, clsize = 10)

  set.seed(5)
  ht <- honest_erf_tree(Y ~ A | X1 + X2 + mod, data = d, modifiers = ~ mod + X1,
                        cluster = ~cl, max_depth = 2, min_samples = 40, n = 21L)

  expect_s3_class(ht, "honest_erf_tree")
  expect_s3_class(ht$split, "factor")
  expect_length(ht$split, nrow(d))
  expect_true(all(tapply(ht$split, d$cl, function(z) length(unique(z))) == 1L))
  expect_identical(ht$tree$size, sum(ht$split == "inference"))
  expect_identical(ht$discovery$tree$size, sum(ht$split == "discovery"))
  expect_true(is.call(ht$call))

  # Without cross-fitting the clustering changes no estimate, only the variance, so the
  # same tree refit on an unclustered ERF of the same rows gives the same node curves with
  # different bands.
  inf <- d[ht$split == "inference", ]
  e_cl <- estimate_erf(Y ~ A | X1 + X2 + mod, data = inf, n = 21L, cluster = ~cl)
  e_iid <- estimate_erf(Y ~ A | X1 + X2 + mod, data = inf, n = 21L)

  rf_cl <- refit_erf_tree(ht$discovery, e_cl)
  rf_iid <- refit_erf_tree(ht$discovery, e_iid)

  est_cl <- rf_cl[[2L]]$erf(rf_cl$values, vcov = TRUE)
  est_iid <- rf_iid[[2L]]$erf(rf_iid$values, vcov = TRUE)

  expect_equal(est_cl$est, est_iid$est)
  expect_false(isTRUE(all.equal(est_cl$vcov, est_iid$vcov)))
})

test_that("honest_erf_tree takes weights by name and forwards the ERF settings", {
  d <- sim_erf_data(n = 600)
  d$w <- runif(nrow(d), .5, 1.5)

  set.seed(6)
  ht <- honest_erf_tree(Y ~ A | X1 + X2 + mod, data = d, modifiers = ~mod,
                        weights = "w", max_depth = 1, min_samples = 40, n = 15L,
                        bw = 1.2, range = .9, metric = "AUC", int_weight_method = "none")

  expect_identical(ht$info$metric, "AUC")
  expect_identical(ht$info$int_weight_method, "none")
  expect_identical(ht$info$max_depth, 1)
  expect_length(ht$values, 15L)
  expect_equal(ht[[1L]]$bw, 1.2)

  # The weights reach the inference ERF, and are those of the inference rows.
  expect_equal(ht$weights, d$w[ht$split == "inference"])

  # A bad tree setting fails before anything is fit.
  expect_error(honest_erf_tree(Y ~ A | X1 + X2 + mod, data = d, metric = "L3"))
  expect_error(honest_erf_tree(Y ~ A | X1 + X2 + mod, data = d, frac = 2))
})

test_that("refit_erf_tree rejects inputs the tree cannot be applied to", {
  d <- sim_erf_data(n = 600)

  set.seed(7)
  h <- fit_halves(d)
  tr <- rpart_erf(h$dis, ~ mod + X1, max_depth = 2, min_samples = 50)

  # A node ERF carries no data to refit on.
  expect_error(refit_erf_tree(tr, tr[[1L]]))

  # An ERF estimated at a different order than the tree was grown from.
  e_d1 <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d[h$split == "inference", ],
                       n = 21L, d = 1, bw = 1.5)
  expect_error(refit_erf_tree(tr, e_d1), "d = 1")

  # An ERF whose data lacks a modifier the tree splits on.
  e_no_mod <- estimate_erf(Y ~ A | X1 + X2, data = d[h$split == "inference", ], n = 21L)
  expect_error(refit_erf_tree(tr, e_no_mod))

  expect_error(refit_erf_tree(h$dis, h$inf))
})

test_that("nodes too small on the inference sample warn or error, by name", {
  d <- sim_erf_data(n = 600)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~mod, max_depth = 1, min_samples = 50)

  expect_identical(tr$tree$split_var, "mod")

  # An inference sample so small that a child falls under min_samples.
  set.seed(8)
  tiny <- d[sample.int(nrow(d), 60L), ]
  e_tiny <- estimate_erf(Y ~ A | X1 + X2 + mod, data = tiny, n = 21L)

  expect_warning(rf <- refit_erf_tree(tr, e_tiny), "fewer than 50 units")
  expect_s3_class(rf, "honest_erf_tree")

  # An inference sample in which one child is empty cannot be refit at all. The modifier
  # is added to the ERF's data after the fit, since a one-level factor cannot enter the
  # nuisance models.
  only0 <- d[d$mod == "0", ]
  e_only0 <- estimate_erf(Y ~ A | X1 + X2, data = only0, n = 21L)
  e_only0$data$mod <- only0$mod

  expect_error(suppressWarnings(refit_erf_tree(tr, e_only0)), "too few")
})

test_that("print, summary, plot, and [[ work on an honest tree", {
  d <- sim_erf_data(n = 600)

  set.seed(9)
  ht <- honest_erf_tree(Y ~ A | X1 + X2 + mod, data = d, modifiers = ~ mod + X1,
                        max_depth = 2, min_samples = 40, n = 21L)

  out <- capture.output(print(ht))
  expect_true(any(grepl("<honest_erf_tree>", out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("discovery sample (n = %d)", ht$discovery$tree$size), out,
                        fixed = TRUE)))
  expect_true(any(grepl(sprintf("inference sample (n = %d)", ht$tree$size), out,
                        fixed = TRUE)))

  # The sizes summary() shows are the inference sizes.
  s <- as.character(summary(ht))
  nodes <- walk_tree(ht$tree)
  expect_length(s, length(nodes))
  for (i in seq_along(nodes)) {
    expect_match(s[i], sprintf("n = %d", nodes[[i]]$size))
  }

  expect_s3_class(ht[[2L]], "erf_node")
  expect_identical(ht[[2L]]$size, ht$tree$left$size)
  expect_s3_class(plot(ht), "ggplot")
  expect_s3_class(plot(ht[[2L]], density = FALSE), "ggplot")

  panels <- plot(ht, type = "list", panels = "marginal")
  expect_length(panels, length(nodes))
})

test_that("rpart_erf node bands are cluster-robust when the ERF is clustered", {
  d <- sim_erf_data(n = 500, clsize = 10)

  e_cl <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L, cluster = ~cl)
  e_iid <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  tr_cl <- rpart_erf(e_cl, ~mod, max_depth = 1, min_samples = 50)
  tr_iid <- rpart_erf(e_iid, ~mod, max_depth = 1, min_samples = 50)

  expect_identical(tr_cl$tree$split_var, tr_iid$tree$split_var)

  est_cl <- tr_cl[[2L]]$erf(tr_cl$values, vcov = TRUE)
  est_iid <- tr_iid[[2L]]$erf(tr_iid$values, vcov = TRUE)

  # Same curve at the node, different bands.
  expect_equal(est_cl$est, est_iid$est)
  expect_false(isTRUE(all.equal(est_cl$vcov, est_iid$vcov)))

  # One unit per cluster collapses the cluster-robust node variance to the iid one.
  d$one_each <- factor(seq_len(nrow(d)))
  e_one <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L, cluster = ~one_each)
  tr_one <- rpart_erf(e_one, ~mod, max_depth = 1, min_samples = 50)

  expect_equal(tr_one[[2L]]$erf(tr_one$values, vcov = TRUE)$vcov, est_iid$vcov)

  # The tree records what it was grown with.
  expect_identical(tr_cl$modifiers, ~mod)
  expect_identical(tr_cl$d, 0L)
})

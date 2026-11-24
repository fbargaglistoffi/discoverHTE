test_that("rpart_erf prefers a real effect modifier to a noise variable", {
  d <- sim_erf_data(n = 500)
  set.seed(99)
  d$noise <- factor(rbinom(nrow(d), 1, 0.5))

  # `noise` has to be in the formula to reach `erf$data`, which is where
  # process_modifiers() resolves candidate modifiers from.
  e <- estimate_erf(Y ~ A | X1 + X2 + mod + noise, data = d, n = 21L)

  tr <- rpart_erf(e, ~ mod + noise, max_depth = 1, min_samples = 50)

  expect_s3_class(tr, "erf_tree")
  expect_identical(tr$tree$type, "internal")
  expect_identical(tr$tree$split_var, "mod")
})

test_that("the modifier is picked up somewhere in a deeper tree", {
  d <- sim_erf_data(n = 500)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  tr <- rpart_erf(e, ~ mod + X1 + X2, max_depth = 3, min_samples = 40)

  split_vars <- function(node) {
    if (node$type == "leaf") return(character(0))
    c(node$split_var, split_vars(node$left), split_vars(node$right))
  }

  expect_true("mod" %in% split_vars(tr$tree))
})

test_that("the tree records the nuisance method and min node size", {
  d <- sim_erf_data()
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  tr <- rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 40)

  expect_identical(attr(tr$method, "method_name"), "glm")
  expect_identical(tr$info$min_samples, 40)

  # Both must actually appear in the printed output; `min_sample` (singular) used to be
  # read here, so the line rendered empty.
  out <- capture.output(print(tr))
  expect_true(any(grepl("Nuisance method: glm", out, fixed = TRUE)))
  expect_true(any(grepl("Min. node size: 40", out, fixed = TRUE)))
})

test_that("max_depth and min_samples bound the tree", {
  d <- sim_erf_data(n = 500)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  leaf_depths <- function(node, depth = 0L) {
    if (node$type == "leaf") return(depth)
    c(leaf_depths(node$left, depth + 1L), leaf_depths(node$right, depth + 1L))
  }
  leaf_sizes <- function(node) {
    if (node$type == "leaf") return(node$size)
    c(leaf_sizes(node$left), leaf_sizes(node$right))
  }

  tr <- rpart_erf(e, ~ mod + X1 + X2, max_depth = 2, min_samples = 60)

  expect_lte(max(leaf_depths(tr$tree)), 2L)
  expect_true(all(leaf_sizes(tr$tree) >= 60))

  # An impossible minimum yields a single leaf.
  stump <- rpart_erf(e, ~ mod + X1, max_depth = 3, min_samples = nrow(d))
  expect_identical(stump$tree$type, "leaf")
  expect_identical(stump$tree$reason, "min_samples")
})

test_that("[[ returns the ERF at a node and node sizes are coherent", {
  d <- sim_erf_data(n = 500)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  tr <- rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 50)

  root <- tr[[1L]]
  expect_s3_class(root, "erf_node")
  expect_s3_class(root, "erf")
  expect_identical(root$size, nrow(d))

  # A parent's size is the sum of its children's.
  expect_identical(tr$tree$size, tr$tree$left$size + tr$tree$right$size)

  expect_error(tr[[999L]], "must be a single node position")
})

test_that("rpart_erf rejects a node ERF, which carries no data", {
  d <- sim_erf_data()
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~ mod + X1, max_depth = 1, min_samples = 50)

  expect_error(rpart_erf(tr[[1L]], ~mod))
})

test_that("split_type picks the candidate cut points", {
  m <- data.frame(v = seq(0, 1, length.out = 1000))

  q <- discoverHTE:::make_splitters("v", m, 5L, "quantile")
  r <- discoverHTE:::make_splitters("v", m, 5L, "range")

  expect_length(q, 4L)
  expect_length(r, 4L)

  # On a uniform grid the two agree; they diverge once the modifier is skewed.
  expect_equal(vapply(q, attr, numeric(1L), "breaks"),
               vapply(r, attr, numeric(1L), "breaks"),
               tolerance = 0.01)

  skewed <- data.frame(v = c(seq(0, 1, length.out = 900), seq(50, 100, length.out = 100)))
  expect_false(identical(
    vapply(discoverHTE:::make_splitters("v", skewed, 5L, "quantile"), attr, numeric(1L), "breaks"),
    vapply(discoverHTE:::make_splitters("v", skewed, 5L, "range"), attr, numeric(1L), "breaks")
  ))
})

test_that("quantile splitting survives a modifier with a concentrated mass", {
  # Duplicated quantiles are dropped rather than passed to cut(), which rejects them.
  set.seed(1)
  ties <- data.frame(v = c(rep(0, 700), runif(300, 1, 10)))

  expect_no_error(sp <- discoverHTE:::make_splitters("v", ties, 5L, "quantile"))
  expect_true(length(sp) >= 1L)

  # When every interior quantile collapses onto one value there is no cut to make, and
  # the node degrades to a leaf rather than erroring.
  allties <- data.frame(v = c(rep(0, 980), runif(20, 1, 10)))
  expect_no_error(sp0 <- discoverHTE:::make_splitters("v", allties, 5L, "quantile"))
  expect_length(sp0, 0L)
})

test_that("split_type is validated and reaches the tree", {
  d <- sim_erf_data(n = 500)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  expect_error(rpart_erf(e, ~ mod + X1, split_type = "deciles"))

  tr_q <- rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 50,
                    split_type = "quantile")
  tr_r <- rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 50,
                    split_type = "range")

  expect_s3_class(tr_q, "erf_tree")
  expect_s3_class(tr_r, "erf_tree")
  expect_identical(tr_q$tree$size, nrow(d))
})

test_that("a modifier with missing values splits rather than erroring", {
  set.seed(3)
  m <- data.frame(v = runif(500))
  m$v[c(7L, 99L, 301L)] <- NA

  for (st in c("quantile", "range")) {
    expect_no_error(sp <- discoverHTE:::make_splitters("v", m, 5L, st))
    expect_true(length(sp) >= 1L)

    # Missing units go right, so the two branches still partition the sample and
    # neither carries an NA that would propagate into the node sizes.
    res <- discoverHTE:::apply_splitter(sp[[1L]], rep(TRUE, nrow(m)), m)
    expect_false(anyNA(res$splitL))
    expect_false(anyNA(res$splitR))
    expect_identical(sum(res$splitL) + sum(res$splitR), nrow(m))
    expect_true(all(res$splitR[is.na(m$v)]))
  }

  # The cut points are read off observed values only, so they are unaffected by the NAs.
  complete <- data.frame(v = m$v[!is.na(m$v)])
  expect_equal(vapply(discoverHTE:::make_splitters("v", m, 5L, "quantile"), attr, numeric(1L), "breaks"),
               vapply(discoverHTE:::make_splitters("v", complete, 5L, "quantile"), attr, numeric(1L), "breaks"))
})

test_that("a tree grows when a modifier has missing values", {
  d <- sim_erf_data(n = 500)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  e$data$X1[c(1L, 2L, 3L)] <- NA

  expect_no_error(tr <- rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 50))
  expect_identical(tr$tree$size, nrow(d))
})

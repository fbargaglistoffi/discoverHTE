# The node positions summary() displays are a documented contract: ?summary.erf_tree and
# ?`[[.erf_tree` both say they come from a depth-first walk counting every node, leaves
# included, and that `tree[[i]]` is the ERF at the node labelled [i]. Nothing tested that
# correspondence, which is easy to break from either side.

# The same walk, written independently of the package's own tree_nodes().
walk_depth_first <- function(node) {
  out <- list(node)

  if (!is.null(node$left)) {
    out <- c(out, walk_depth_first(node$left))
  }

  if (!is.null(node$right)) {
    out <- c(out, walk_depth_first(node$right))
  }

  out
}

test_that("summary() numbers nodes depth-first, counting leaves", {
  d <- sim_erf_data(n = 600)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~ mod + X1, max_depth = 3, min_samples = 40)

  nodes <- walk_depth_first(tr$tree)
  out <- as.character(summary(tr))

  expect_length(out, length(nodes))

  # Each displayed line carries its position and that node's size, in walk order.
  for (i in seq_along(nodes)) {
    expect_match(out[i], sprintf("\\[%d\\]", i))
    expect_match(out[i], sprintf("n = %d", nodes[[i]]$size), fixed = FALSE)
  }
})

test_that("summary() marks exactly the leaves", {
  d <- sim_erf_data(n = 600)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~ mod + X1, max_depth = 3, min_samples = 40)

  nodes <- walk_depth_first(tr$tree)
  out <- as.character(summary(tr))

  starred <- grepl("\\*\\s*$", out)
  is_leaf <- vapply(nodes, function(n) identical(n$type, "leaf"), logical(1L))

  expect_identical(starred, is_leaf)
  expect_true(any(is_leaf))
  expect_false(all(is_leaf))
})

test_that("`[[` returns the node at the position summary() displays", {
  d <- sim_erf_data(n = 600)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~ mod + X1, max_depth = 3, min_samples = 40)

  nodes <- walk_depth_first(tr$tree)

  for (i in seq_along(nodes)) {
    node_erf <- tr[[i]]

    expect_s3_class(node_erf, "erf_node")
    expect_identical(node_erf$size, nodes[[i]]$size)
    expect_identical(node_erf$path, nodes[[i]]$path)
  }

  # A string still indexes the list, as for any other list.
  expect_identical(tr[["values"]], tr$values)
  expect_type(tr[["tree"]], "list")
})

test_that("`[[` rejects a position outside the tree", {
  d <- sim_erf_data(n = 400)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 50)

  n_nodes <- length(walk_depth_first(tr$tree))

  expect_error(tr[[0L]], "single node position")
  expect_error(tr[[n_nodes + 1L]], "single node position")
  expect_error(tr[[c(1L, 2L)]], "single node position")
  expect_error(tr[[NA_integer_]], "single node position")

  # The last valid position works, so the bound is not off by one.
  expect_s3_class(tr[[n_nodes]], "erf_node")
})

test_that("a stump has exactly one node at position 1", {
  d <- sim_erf_data(n = 300)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)

  stump <- rpart_erf(e, ~ mod + X1, max_depth = 3, min_samples = nrow(d))

  expect_identical(stump$tree$type, "leaf")
  expect_length(as.character(summary(stump)), 1L)
  expect_identical(stump[[1L]]$size, nrow(d))
  expect_error(stump[[2L]], "single node position")
})

test_that("print() on a tree reports its settings and summary() returns its own class", {
  d <- sim_erf_data(n = 400)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~ mod + X1, metric = "L1", int_weight_method = "gps",
                  max_depth = 2, min_samples = 50)

  s <- summary(tr)
  expect_s3_class(s, "summary.erf_tree")
  expect_type(unclass(s), "character")

  out <- capture.output(print(tr))
  # The label follows the `modifiers` argument; it read "Moderators:" before.
  expect_true(any(grepl("Modifiers: mod and X1", out, fixed = TRUE)))
  expect_false(any(grepl("Moderators", out, fixed = TRUE)))
  expect_true(any(grepl("Splitting metric: L1", out, fixed = TRUE)))
  expect_true(any(grepl("Metric weight: gps", out, fixed = TRUE)))
  expect_true(any(grepl("Max. tree depth: 2", out, fixed = TRUE)))

  # print() returns its argument invisibly.
  capture.output(vis <- withVisible(print(tr))$visible)
  expect_false(vis)
})

test_that("a node ERF prints as a node, not as a full erf", {
  d <- sim_erf_data(n = 400)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  tr <- rpart_erf(e, ~ mod + X1, max_depth = 2, min_samples = 50)

  out <- capture.output(print(tr[[2L]]))

  expect_true(any(grepl("one node", out, fixed = TRUE)))
  expect_true(any(grepl("Node size:", out, fixed = TRUE)))

  # The root has no branch path and says so.
  expect_true(any(grepl("root", capture.output(print(tr[[1L]])), fixed = TRUE)))
})

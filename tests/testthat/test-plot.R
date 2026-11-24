# ?plot.erf_tree documents several contracts that nothing tested: which nodes get a panel
# under each `panels` value, that a list's names are the node positions summary() shows,
# and that a marginal panel draws the same curve as plot(x[[i]]).

walk_nodes <- function(node) {
  out <- list(node)

  if (!is.null(node$left)) {
    out <- c(out, walk_nodes(node$left))
  }

  if (!is.null(node$right)) {
    out <- c(out, walk_nodes(node$right))
  }

  out
}

fit_tree <- function(n = 600, max_depth = 3, min_samples = 40) {
  d <- sim_erf_data(n = n)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  rpart_erf(e, ~ mod + X1, max_depth = max_depth, min_samples = min_samples)
}

test_that("plot.erf returns a ggplot for every band and density setting", {
  d <- sim_erf_data(n = 400)
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)

  expect_s3_class(plot(e, density = FALSE, simultaneous = FALSE), "ggplot")
  expect_s3_class(plot(e, conf_level = 0, density = FALSE), "ggplot")

  skip_if_not_installed("mvtnorm")
  expect_s3_class(plot(e, density = FALSE, simultaneous = TRUE), "ggplot")

  skip_if_not_installed("ggside")
  expect_s3_class(plot(e), "ggplot")
})

test_that("plot.erf drops the ribbon when no bands are asked for", {
  d <- sim_erf_data(n = 400)
  e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 21L)

  with_band <- plot(e, density = FALSE, simultaneous = FALSE)
  no_band <- plot(e, conf_level = 0, density = FALSE)

  expect_length(with_band$layers, 2L)
  expect_length(no_band$layers, 1L)

  expect_error(plot(e, conf_level = 1, density = FALSE))
  expect_error(plot(e, conf_level = -0.5, density = FALSE))
})

test_that("type = 'tree' gives one plot and type = 'list' gives named panels", {
  tr <- fit_tree()

  whole <- plot(tr, type = "tree")
  expect_s3_class(whole, "ggplot")

  panels <- plot(tr, type = "list")
  expect_type(panels, "list")
  expect_true(all(vapply(panels, inherits, logical(1L), "ggplot")))

  # "nodes" is kept as a synonym for "list".
  expect_identical(names(plot(tr, type = "nodes")), names(panels))
})

test_that("panels = 'children' panels exactly the internal nodes", {
  tr <- fit_tree()

  nodes <- walk_nodes(tr$tree)
  internal <- which(vapply(nodes, function(n) identical(n$type, "internal"), logical(1L)))

  panels <- plot(tr, type = "list", panels = "children")

  expect_identical(names(panels), as.character(internal))
})

test_that("panels = 'marginal' panels every node, leaves included", {
  tr <- fit_tree()

  nodes <- walk_nodes(tr$tree)
  panels <- plot(tr, type = "list", panels = "marginal")

  expect_identical(names(panels), as.character(seq_along(nodes)))
})

test_that("each panel carries the path of the node it belongs to", {
  tr <- fit_tree()

  nodes <- walk_nodes(tr$tree)
  panels <- plot(tr, type = "list", panels = "marginal")

  for (i in seq_along(nodes)) {
    expect_identical(attr(panels[[i]], "path"), nodes[[i]]$path)
  }
})

test_that("a marginal panel draws the same curve as plot(x[[i]])", {
  # Documented in ?`[[.erf_tree`: plot(x, type = "list", panels = "marginal") gives the
  # same per-node ERFs as x[[i]].
  tr <- fit_tree()

  panels <- plot(tr, type = "list", panels = "marginal")

  for (i in c(1L, 2L, length(panels))) {
    from_panel <- ggplot2::layer_data(panels[[i]], 1L)
    from_node <- ggplot2::layer_data(plot(tr[[i]], conf_level = 0, density = FALSE), 1L)

    expect_equal(from_panel$y, from_node$y, tolerance = 1e-8)
  }
})

test_that("plot.erf_node subtitles a node with its branch path, except the root", {
  tr <- fit_tree(max_depth = 2, min_samples = 50)

  root <- plot(tr[[1L]], conf_level = 0, density = FALSE)
  expect_null(root$labels$subtitle)

  child <- plot(tr[[2L]], conf_level = 0, density = FALSE)
  expect_type(child$labels$subtitle, "character")

  # The subtitle is the node's own branch path, joined with `&`.
  expect_identical(child$labels$subtitle,
                   paste(tr$tree$left$path[-1L], collapse = " & "))
})

test_that("plot.erf_tree rejects arguments outside its allowed values", {
  tr <- fit_tree(max_depth = 2, min_samples = 50)

  expect_error(plot(tr, type = "sideways"))
  expect_error(plot(tr, panels = "siblings"))
})

test_that("a stump still plots", {
  d <- sim_erf_data(n = 300)
  e <- estimate_erf(Y ~ A | X1 + X2 + mod, data = d, n = 21L)
  stump <- rpart_erf(e, ~ mod + X1, max_depth = 3, min_samples = nrow(d))

  # With no internal node there is nothing to panel under "children", but the marginal
  # view still has the root.
  marginal <- plot(stump, type = "list", panels = "marginal")
  expect_length(marginal, 1L)
  expect_s3_class(marginal[[1L]], "ggplot")
})

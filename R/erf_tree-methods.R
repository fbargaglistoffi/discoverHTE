#' Extract from an ERF Tree
#'
#' @description
#' `[[` extracts a component of an `<erf_tree>` by name, as for any list, but extracts
#' the **ERF at a node** when given a node position. The positions are those shown by
#' [summary()][summary.erf_tree] (a depth-first walk of the tree, counting every node,
#' leaves included), so `tree[[2]]` is the ERF of the node labelled `[2]` there.
#'
#' @param x an `<erf_tree>` object; the output of a call to [rpart_erf()].
#' @param i a string naming a component of `x`, or a single number giving a node
#'   position.
#' @param ... ignored.
#'
#' @return
#' For a string, the named component. For a number, a thin `<erf>` object for that node,
#' which can be plotted with [plot.erf()].
#'
#' @details
#' The ERF stored at a node is a *thin* `<erf>`: it holds the node's estimates and
#' covariance at the anchor values, the exposure values of the units in that node, and
#' the bandwidth. [plot.erf()] needs no more than that, and storing no more avoids carrying
#' a copy of the data and the fitted nuisance models at every node. Evaluating it returns
#' exactly what the full ERF would; [update()][update.erf] does not work on it, though,
#' because re-smoothing requires the nuisance fit.
#'
#' Note the distinction from [plot.erf_tree()]: `plot(x[[i]])` draws the single ERF *at*
#' node `i`, whereas `plot(x, type = "list")[[i]]` draws the two ERFs of the *children* of
#' the `i`-th internal node. `plot(x, type = "list", panels = "marginal")` gives the same
#' per-node ERFs as `x[[i]]`, as panels sharing one y-axis.
#'
#' @seealso [rpart_erf()], [plot.erf()], [summary.erf_tree()]
#'
#' @exportS3Method `[[` erf_tree
`[[.erf_tree` <- function(x, i, ...) {
  if (!is.numeric(i)) {
    return(NextMethod())
  }

  nodes <- tree_nodes(x$tree)

  if (length(i) != 1L || anyNA(i) || i < 1 || i > length(nodes)) {
    arg::err("{.arg i} must be a single node position between 1 and {length(nodes)}; use {.fn summary} to see the tree")
  }

  nodes[[i]][["erf"]]
}

# The nodes of a tree in the same depth-first order as node_ids(), i.e. the order in
# which extract_tree() and summary() number them.
tree_nodes <- function(x) {
  out <- list(x)

  if (!is_null(x$left)) {
    out <- c(out, Recall(x$left))
  }

  if (!is_null(x$right)) {
    out <- c(out, Recall(x$right))
  }

  out
}

#' Plot the ERF at a Node of an ERF Tree
#'
#' @description
#' Plots the ERF held at one node of an `<erf_tree>`, as extracted by
#' [`[[.erf_tree`]. This is [plot.erf()] with the node's branch path added as a
#' subtitle, so a node ERF plotted on its own still says which subgroup it describes.
#'
#' @param x an `<erf>` object for a single node, e.g. `tree[[2]]`.
#' @param ... passed to [plot.erf()].
#'
#' @return
#' A `<ggplot>` object.
#'
#' @seealso [plot.erf()], [`[[.erf_tree`]
#'
#' @exportS3Method plot erf_node
plot.erf_node <- function(x, ...) {
  p <- NextMethod()

  subtitle <- path_label(x[["path"]])

  if (is_null(subtitle)) {
    return(p)
  }

  p + labs(subtitle = subtitle)
}

#' @exportS3Method print erf_node
print.erf_node <- function(x, ...) {
  cli::format_inline("An {.cls erf} object for one node of an {.cls erf_tree}\n") |>
    cli::style_italic() |>
    cli::cat_line()

  node <- path_label(x[["path"]]) %or% "root"

  treat_name <- attr(x$treat, "varname") %or% "exposure"

  cli::format_inline(" - Node: {node}\n",
                     " - Node size: {x$size}\n",
                     " - Exposure: {treat_name}\n",
                     " - Exposure range: {round(min(x$values), 4)} to {round(max(x$values), 4)}\n",
                     " - Bandwidth: {round(x$bw, 4)}") |>
    cli::cat_line()

  invisible(x)
}

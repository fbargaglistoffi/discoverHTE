#' Summarize and Display an ERF Tree
#' @name summary.erf_tree
#'
#' @description
#' `print()` displays the settings used to grow an `<erf_tree>`, and `summary()` displays
#' the tree itself, with every node numbered and labeled by the split it descends from.
#' For an `<honest_erf_tree>`, `print()` also gives the sizes of the discovery and
#' inference samples, and the node sizes `summary()` shows are those on the inference
#' sample.
#'
#' @param x,object an `<erf_tree>` object; the output of a call to [rpart_erf()] or
#'   [honest_erf_tree()]. For `print()` on the output of `summary()`, a
#'   `<summary.erf_tree>` object.
#' @param ... ignored.
#'
#' @return
#' For `summary()`, a `<summary.erf_tree>` object, which is a character vector with one
#' element per line of the displayed tree and a `print()` method that writes those lines to
#' the console. For `print()`, its argument, invisibly.
#'
#' @details
#' `summary()` draws the tree as nested branches. Each node is prefixed with its position in
#' square brackets and followed by its size, and a leaf is marked with `*`. The positions
#' come from a depth-first walk that counts every node, leaves included, so they descend the
#' left branch of a node before crossing to the right. These are the positions
#' [`[[.erf_tree`] takes to extract the ERF at a node, and the names given to the panels of
#' `plot(x, type = "list")`.
#'
#' Every line below the root is labeled with the split the node descends from (e.g.,
#' `PIR < 2.77`), so the subgroup a node describes is read by following the labels from the
#' root down to it.
#'
#' @seealso [rpart_erf()], [`[[.erf_tree`], [plot.erf_tree()]
#'
#' @examples
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
#' tr <- rpart_erf(e, ~ M + X1, max_depth = 2, min_samples = 50)
#'
#' # The settings the tree was grown with
#' tr
#'
#' # The tree itself, each node numbered and labeled by its split
#' summary(tr)
#'
#' # Those numbers are the positions `[[` takes
#' tr[[2]]
#'
#' @exportS3Method print erf_tree
print.erf_tree <- function(x, ...) {
  # out <- extract_tree(x$tree)

  cli::format_inline("An {.cls {class(x)[1L]}} object\n") |>
    cli::style_italic() |>
    cli::cat_line()

  nm <- names(x$M)
  method <- attr(x$method, "method_name") %or% "unknown"

  cli::format_inline(" - Modifiers: {nm}\n",
                     " - Nuisance method: {method}\n",
                     " - Splitting metric: {x$info$metric}\n",
                     " - Metric weight: {x$info$int_weight_method}\n",
                     " - Max. tree depth: {x$info$max_depth}\n",
                     " - Min. node size: {x$info$min_samples}") |>
    cli::cat_line()

  # An honest tree was grown on one sample and had its node ERFs estimated on another;
  # the root sizes give the two sample sizes.
  if (inherits(x, "honest_erf_tree")) {
    cli::format_inline(" - Grown on: discovery sample (n = {x$discovery$tree$size})\n",
                       " - Node ERFs: inference sample (n = {x$tree$size})") |>
      cli::cat_line()
  }

  cli::cat_line()

  if (interactive()) {
    width <- min(60, getOption("width", 60))
    sprintf("{.emph Use {.help [{.fun summary}](%1$s::summary.erf_tree)} to display the tree.}",
            utils::packageName()) |>
      cli::format_inline() |>
      cli::ansi_strwrap(width) |>
      cli::cat_line()
  }

  invisible(x)
}

#' @rdname summary.erf_tree
#' @exportS3Method summary erf_tree
summary.erf_tree <- function(object, ...) {
  out <- extract_tree(object$tree)
  class(out) <- "summary.erf_tree"
  out
}

#' @rdname summary.erf_tree
#' @exportS3Method print summary.erf_tree
print.summary.erf_tree <- function(x, ...) {
  x |> unclass() |>
    cli::cat_line()

  invisible(x)
}

node_ids <- function(x) {
  out <- list(x[["path"]])

  if (!is_null(x$left)) {
    out <- c(out, Recall(x$left))
  }

  if (!is_null(x$right)) {
    out <- c(out, Recall(x$right))
  }

  out
}

extract_tree <- function(x, ids = NULL, bars = logical(0), first = TRUE) {
  if (is_null(ids)) {
    ids <- node_ids(x)
  }

  names <- as.character(seq_along(ids))

  tree_path <- c()

  if (first) {
    tree_path <- sprintf("[%s] root (n = %s)",
                         names[which(vapply(ids, identical, logical(1L), x$path))],
                         x$size)

    if (x$type == "leaf") {
      tree_path <- paste(tree_path, "*")
    }
  }

  if (x$type != "leaf") {

    slabs <- x$split_levels |>
      setNames(c("left", "right"))

    prefix <- {
      if (is_null(bars)) " "
      else paste0(" ", paste(ifelse(bars, "\u2502    ", "     "), collapse = ""))
    }

    for (i in c("left", "right")) {
      node <- x[[i]]

      knam <- names[which(vapply(ids, identical, logical(1L), node$path))]

      klab <- if (node$type == "leaf") " *" else ""

      lab <- sprintf("%s%s\u2500\u2500 [%s] %s (n = %s)%s",
                     prefix, switch(i, "right" = "\u2514", "\u251c"), knam, slabs[i], node$size, klab)

      bars_new <- c(bars, i != "right")

      tree_path <- c(tree_path,
                     lab,
                     extract_tree(node, bars = bars_new, ids = ids,
                                  first = FALSE))
    }
  }

  tree_path
}



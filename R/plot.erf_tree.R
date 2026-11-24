#' Plot an ERF Tree
#'
#' @description
#' Plots the subgroup dose-response functions discovered by [rpart_erf()]. Two choices
#' are independent of one another: what each panel shows (`panels`), and whether the
#' panels are laid out as a tree or returned as a list (`type`).
#'
#' @param x an `<erf_tree>` object; the output of a call to [rpart_erf()].
#' @param type the shape of the output. `"tree"` (the default) arranges the panels in a
#'   single plot laid out as the tree, with labelled branches connecting a node to its
#'   children. `"list"` returns a list of the same panels, named by the node positions
#'   shown by [summary()][summary.erf_tree]. `"nodes"` is accepted as a synonym for
#'   `"list"`.
#' @param panels what each panel shows. `"children"` (the default) draws one ERF per
#'   child, so the two lines in a panel are the dose-responses either side of that
#'   node's split; only internal nodes get a panel. `"marginal"` draws the single ERF at
#'   the node itself, so *every* node gets a panel, leaves included. A marginal tree is
#'   therefore wider and more crowded than a children tree of the same depth.
#' @param conf_level the confidence level for the pointwise bands. Default is `.95`; set
#'   to `0` to omit the bands.
#' @param fixed_y whether to use a common y-axis across panels. Default is `TRUE`, which
#'   makes the panels directly comparable.
#' @param panel_scale for `type = "tree"`, the fraction of the available width each
#'   panel fills. Default is `0.9`; lower values leave more room between panels.
#' @param branch_space for `type = "tree"`, the fraction of each row left free for the
#'   branches and their labels. Default is `0.22`.
#' @param label_size for `type = "tree"`, the text size of the branch labels.
#' @param label_width for `type = "tree"`, the width in characters at which branch
#'   labels are wrapped. Default is `24`.
#' @param ... ignored.
#'
#' @return
#' For `type = "tree"`, a single `ggplot` object. For `type = "list"`, a list of `ggplot`
#' objects, each with a `"path"` attribute identifying its node, named by the node
#' positions shown by [summary()][summary.erf_tree]. With `panels = "children"` only
#' internal nodes appear, so those names (rather than the list positions) line up with
#' [`[[.erf_tree`].
#'
#' @details
#' With `panels = "children"`, a terminal leaf gets no panel of its own: its
#' dose-response is one of the two lines in its parent's panel. Branches are therefore
#' drawn only to children that are themselves internal nodes, and the branch to a leaf
#' stops at its label. With `panels = "marginal"` every node has a panel, so every branch
#' ends in one.
#'
#' Panels are positioned according to their place in the tree, with the layout derived
#' from the tree's own shape, so it adapts to any structure.
#'
#' `plot(x, type = "list", panels = "marginal")[[i]]` and `plot(x[[i]])` both draw the
#' ERF at node `i`. The former uses this function's pointwise bands and shared y-axis so
#' that it sits alongside the other panels; the latter is the fuller standalone view from
#' [plot.erf()], with simultaneous bands and an exposure density by default.
#'
#' @seealso [rpart_erf()], [plot.erf()], [`[[.erf_tree`]
#'
#' @exportS3Method plot erf_tree
plot.erf_tree <- function(x, type = "tree", panels = "children", conf_level = .95,
                          fixed_y = TRUE, panel_scale = 0.9, branch_space = 0.22,
                          label_size = 2.8, label_width = 24L, ...) {
  type <- arg::match_arg(type, c("tree", "list", "nodes"))
  panels <- arg::match_arg(panels, c("children", "marginal"))

  arg::arg_number(conf_level)
  arg::arg_between(conf_level, c(0, 1), inclusive = c(TRUE, FALSE))

  if (type == "tree") {
    return(plot_erf_tree_native(x, panels = panels, conf_level = conf_level,
                                fixed_y = fixed_y, panel_scale = panel_scale,
                                branch_space = branch_space, label_size = label_size,
                                label_width = label_width))
  }

  # The nodes in the order summary() numbers them, so the list can be named by those
  # positions.
  nodes <- tree_nodes(x$tree)

  keep <- {
    if (panels == "marginal") seq_along(nodes)
    else which(vapply(nodes, function(nd) identical(nd$type, "internal"), logical(1L)))
  }

  if (is_null(keep)) {
    arg::err("the tree has no internal nodes, so there is nothing to plot; use {.code panels = \"marginal\"} to plot the ERF at the root")
  }

  node_data <- lapply(nodes[keep], panel_data, values = x$values,
                      conf_level = conf_level, panels = panels)

  ylims <- {
    if (!fixed_y) NULL
    else frange(unlist(lapply(node_data, function(d) c(d$est, d$ci_lo, d$ci_hi))))
  }

  treat_name <- get_treat_name(x)

  plots <- lapply(seq_along(keep), function(i) {
    path <- nodes[[keep[i]]][["path"]]

    p <- make_node_panel(node_data[[i]], treat_name = treat_name,
                         conf_level = conf_level, ylims = ylims,
                         subtitle = path_label(path))

    attr(p, "path") <- path

    p
  })

  names(plots) <- as.character(keep)

  plots
}

# The exposure name, taken from the tree object when available. Note that `x$tree$erf` is
# the node's ERF *function* for trees built by rpart_erf(), so it cannot be subset with
# `$`; only fall back to it when it is a list carrying `treat` (as produced by the
# development honest-tree builders).
get_treat_name <- function(x) {
  x[["treat_name"]] %or%
    (if (!is.function(x$tree$erf)) attr(x$tree$erf$treat, "varname")) %or%
    "exposure"
}

# One node's panel: an ERF line (and band) per child, labelled by the split levels.
# `compact` shrinks the text and legend so the panel stays readable when it is one cell
# of a tree layout.
make_node_panel <- function(d, treat_name, conf_level = .95, ylims = NULL,
                            compact = FALSE, subtitle = NULL) {
  # A children panel carries one curve per child, told apart by `group`; a marginal panel
  # carries the single ERF at the node and so needs no colour and no legend.
  grouped <- "group" %in% names(d)

  p <- ggplot(data = d) +
    {
      if (grouped) geom_line(aes(y = .data$est, x = .data$treat, color = .data$group))
      else geom_line(aes(y = .data$est, x = .data$treat))
    }

  if (conf_level > 0) {
    p <- p +
      {
        if (grouped) geom_ribbon(aes(ymin = .data$ci_lo, ymax = .data$ci_hi,
                                     x = .data$treat, fill = .data$group), alpha = .2)
        else geom_ribbon(aes(ymin = .data$ci_lo, ymax = .data$ci_hi,
                             x = .data$treat), alpha = .2)
      }
  }

  # NB: `expand` must be a single logical per side; a bare `c(TRUE, FALSE)` errors on
  # ggplot2 >= 4.x ("if (!expand) ... condition has length > 1").
  p <- p +
    coord_cartesian(expand = c(left = TRUE, right = TRUE), default = TRUE,
                    ylim = ylims)

  p <- p +
    theme_bw(base_size = if (compact) 8 else 11) +
    theme(legend.position = if (grouped) "bottom" else "none",
          legend.title = element_blank()) +
    labs(x = treat_name,
         y = "E[Y]",
         # The tree layout already shows each panel's path, and its panels are small.
         subtitle = if (!compact) subtitle)

  if (compact) {
    p <- p +
      theme(legend.margin = margin(0, 0, 0, 0),
            legend.box.spacing = unit(2, "pt"),
            legend.key.size = unit(7, "pt"),
            legend.text = element_text(size = 6),
            plot.margin = margin(1, 3, 1, 1))
  }

  p
}

get_erf_est <- function(erf, values, conf_level = .95) {
  est <- erf(values, vcov = conf_level > 0)

  data.frame(
    treat = values,
    est = est$est,
    se = if (conf_level > 0) sqrt(diag(est$vcov)) else 0
  )
}

# The data behind one panel: the ERFs of a node's two children, or the node's own ERF.
panel_data <- function(node, values, conf_level = .95, panels = "children") {
  if (panels == "marginal") {
    get_marginal_est(node, values, conf_level)
  }
  else {
    get_node_est(node, values, conf_level)
  }
}

# The node's own ERF, with the same pointwise band convention as get_node_est().
get_marginal_est <- function(node, values, conf_level = .95) {
  d <- get_erf_est(node_erf_fun(node), values, conf_level)

  if (conf_level > 0) {
    d$ci_lo <- d$est - sqrt(qchisq(conf_level, 1)) * d$se
    d$ci_hi <- d$est + sqrt(qchisq(conf_level, 1)) * d$se
  }

  d
}

get_node_est <- function(node, values, conf_level = .95) {
  d <- rbind(
    get_erf_est(node_erf_fun(node$left), values, conf_level) |>
      cbind(group = 1L),
    get_erf_est(node_erf_fun(node$right), values, conf_level) |>
      cbind(group = 2L)
  )

  d$group <- factor(d$group, levels = 1:2,
                    labels = node$split_levels)

  if (conf_level > 0) {
    d$ci_lo <- d$est - sqrt(qchisq(conf_level, 1)) * d$se
    d$ci_hi <- d$est + sqrt(qchisq(conf_level, 1)) * d$se
  }

  d
}

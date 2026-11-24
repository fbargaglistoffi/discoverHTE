# Native tree layout for <erf_tree> plots -------------------------------------
#
# Panels are produced by make_node_panel() and placed at the position matching their
# location in the tree, computed from the tree's actual shape and composed onto a blank
# ggplot with annotation_custom(), so nothing beyond ggplot2 is needed.
#
# Which nodes get a panel depends on `panels`. With "children" a panel shows one ERF line
# per child, so only internal nodes have one and the panels form a tree over the internal
# nodes: a branch runs to a child's panel only when that child is itself internal, and
# the branch to a leaf stops at its label. With "marginal" a panel shows the node's own
# ERF, so every node has one and every branch ends in a panel.

# Flatten the tree to the nodes that get a panel, in pre-order, recording each one's
# parent and the branch of the parent it descends from. With `panels = "children"` only
# internal nodes get a panel, because a leaf has no children to compare; with
# `panels = "marginal"` every node does.
collect_panel_nodes <- function(tree, panels = "children") {
  out <- list()

  walk <- function(node, parent_id, side) {
    if (is_null(node) ||
        (panels != "marginal" && !identical(node[["type"]], "internal"))) {
      return(invisible(NULL))
    }

    id <- length(out) + 1L

    out[[id]] <<- list(node = node,
                       parent = parent_id,
                       side = side)

    walk(node[["left"]], id, 1L)
    walk(node[["right"]], id, 2L)

    invisible(NULL)
  }

  walk(tree, NA_integer_, NA_integer_)

  out
}

# Tidy layout over the internal-node tree. Panels with no internal children take
# successive slots left to right; every other panel is centred over its children. Rows
# are evenly spaced, deepest at the bottom.
#
# When a node has exactly one internal child, its other (leaf) child is given a position
# of its own as a "stub": the branch and its label are drawn, with no panel at the end.
# That keeps the parent centred between its two branches and keeps the label visible for
# the side that stops at a leaf. A node whose children are *both* leaves gets no stubs --
# it has nothing below it, and its own panel already labels both branches.
erf_tree_layout <- function(nodes) {
  n <- length(nodes)

  parent <- vapply(nodes, function(z) as.integer(z[["parent"]]), integer(1L))
  side <- vapply(nodes, function(z) as.integer(z[["side"]]), integer(1L))

  level <- integer(n)
  for (i in seq_len(n)) {
    if (!is.na(parent[i])) {
      level[i] <- level[parent[i]] + 1L
    }
  }

  kids <- lapply(seq_len(n), function(i) which(!is.na(parent) & parent == i))

  x <- rep(NA_real_, n)
  stub_parent <- integer(0)
  stub_side <- integer(0)
  stub_slot <- numeric(0)
  slot <- 0

  place <- function(i) {
    ks <- kids[[i]]

    if (is_null(ks)) {
      slot <<- slot + 1
      x[i] <<- slot
      return(invisible(NULL))
    }

    child_x <- numeric(0)

    # Left branch then right branch, so slots increase from left to right.
    for (sd in 1:2) {
      k <- ks[side[ks] == sd]

      if (!is_null(k)) {
        place(k)
        child_x <- c(child_x, x[k])
      }
      else if (length(ks) == 1L) {
        slot <<- slot + 1
        stub_parent <<- c(stub_parent, i)
        stub_side <<- c(stub_side, sd)
        stub_slot <<- c(stub_slot, slot)
        child_x <- c(child_x, slot)
      }
    }

    x[i] <<- mean(child_x)

    invisible(NULL)
  }

  place(1L)

  n_rows <- max(level) + 1L

  row_y <- function(l) 1 - (l + 0.5) / n_rows

  list(panels = data.frame(id = seq_len(n),
                           parent = parent,
                           side = side,
                           level = level,
                           # slot centres, rescaled to (0, 1)
                           x = (x - 0.5) / slot,
                           y = row_y(level)),
       stubs = data.frame(parent = stub_parent,
                          side = stub_side,
                          x = (stub_slot - 0.5) / slot,
                          y = row_y(level[stub_parent] + 1L)),
       n_cols = slot,
       n_rows = n_rows)
}

# Compose the panels, branches, and branch labels into a single ggplot.
plot_erf_tree_native <- function(x, panels = "children", conf_level = .95,
                                 fixed_y = TRUE, panel_scale = 0.9, branch_space = 0.22,
                                 label_size = 2.8, label_width = 24L, ...) {
  nodes <- collect_panel_nodes(x[["tree"]], panels = panels)

  if (is_null(nodes)) {
    arg::err("the tree has no internal nodes, so there is nothing to plot; use {.code panels = \"marginal\"} to plot the ERF at the root")
  }

  arg::arg_number(panel_scale)
  arg::arg_between(panel_scale, c(0, 1), inclusive = c(FALSE, TRUE))
  arg::arg_number(branch_space)
  arg::arg_between(branch_space, c(0, 1), inclusive = c(TRUE, FALSE))

  lay <- erf_tree_layout(nodes)
  pan <- lay[["panels"]]
  stubs <- lay[["stubs"]]

  treat_name <- get_treat_name(x)

  node_data <- lapply(nodes, function(z) {
    panel_data(z[["node"]], x[["values"]], conf_level, panels = panels)
  })

  ylims <- {
    if (!fixed_y) NULL
    else frange(unlist(lapply(node_data, function(d) c(d$est, d$ci_lo, d$ci_hi))))
  }

  # Panel extent, in the same units as the layout coordinates. Each row keeps
  # `branch_space` of its height free so the branches and their labels have somewhere to
  # go; without it the panels of adjacent rows very nearly touch.
  w <- panel_scale / lay[["n_cols"]]
  h <- (1 - branch_space) / lay[["n_rows"]]

  # Branches run from the bottom of a parent's panel to the top of where the child sits,
  # and carry the label of the parent's branch that the child descends from. Stub
  # branches are drawn the same way; they simply have no panel at the far end.
  branch_label <- function(parent_id, sd) {
    nodes[[parent_id]][["node"]][["split_levels"]][sd]
  }

  branch <- function(parent_id, side, xend, yend, stub) {
    data.frame(x = pan$x[parent_id],
               y = pan$y[parent_id] - h / 2,
               xend = xend,
               yend = yend + h / 2,
               side = side,
               stub = stub,
               label = vapply(seq_along(parent_id), function(j) {
                 branch_label(parent_id[j], side[j])
               }, character(1L)))
  }

  kid_rows <- which(!is.na(pan$parent))

  edges <- rbind(
    if (is_null(kid_rows)) NULL
    else branch(pan$parent[kid_rows], pan$side[kid_rows], pan$x[kid_rows], pan$y[kid_rows],
                stub = FALSE),
    if (nrow(stubs) == 0L) NULL
    else branch(stubs$parent, stubs$side, stubs$x, stubs$y, stub = TRUE)
  )

  if (!is_null(edges)) {
    # Labels sit at the midpoint of their branch, so a node's two labels line up. Long
    # labels (e.g. a list of factor levels) are wrapped, which is what keeps two wide
    # sibling labels from colliding.
    edges$xlab <- (edges$x + edges$xend) / 2
    edges$ylab <- (edges$y + edges$yend) / 2

    edges$label <- vapply(edges$label, function(l) {
      paste(strwrap(l, width = label_width), collapse = "\n")
    }, character(1L), USE.NAMES = FALSE)

    # A stub branch leads to a leaf, so it stops at its label instead of carrying on to
    # where a panel would have been.
    edges$xend[edges$stub] <- edges$xlab[edges$stub]
    edges$yend[edges$stub] <- edges$ylab[edges$stub]
  }

  out <- ggplot()

  if (!is_null(edges)) {
    out <- out +
      geom_segment(data = edges,
                   aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                   colour = "grey40")
  }

  # annotation_custom() places each panel at its own position; the grobs do not train
  # the scales, so the panel limits are set explicitly below.
  for (i in seq_len(nrow(pan))) {
    panel <- make_node_panel(node_data[[i]], treat_name = treat_name,
                             conf_level = conf_level, ylims = ylims,
                             compact = TRUE)

    out <- out +
      annotation_custom(ggplotGrob(panel),
                        xmin = pan$x[i] - w / 2, xmax = pan$x[i] + w / 2,
                        ymin = pan$y[i] - h / 2, ymax = pan$y[i] + h / 2)
  }

  if (!is_null(edges)) {
    out <- out +
      geom_label(data = edges,
                 aes(x = .data$xlab, y = .data$ylab, label = .data$label),
                 size = label_size, fill = "white", linewidth = 0.25,
                 lineheight = 0.9)
  }

  out +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
    theme_void()
}

get_varnames <- function(expr) {
  recurse <- function(e) {
    if (is.symbol(e)) {
      # bare variable like age
      return(as.character(e))
    }

    if (!is.call(e)) {
      return(NULL)
    }

    # keep as-is for $, [[, and [
    fn <- e[[1L]]
    if (fn == as.name("$") || fn == as.name("[[") || fn == as.name("[")) {
      return(deparse1(e))
    }

    # strip outer function, recurse into arguments
    lapply(as.list(e)[-1L], recurse) |>
      unlist()
  }

  recurse(expr)
}

is_null <- function(x) {isTRUE(length(x) == 0L)}
null_or_error <- function(x) {is_null(x) || inherits(x, "try-error")}
`%or%` <- function(x, y) {
  if (is_null(x)) y else x
}
check_if_zero <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  # this is the default tolerance used in all.equal
  abs(x) < tolerance
}

round_midpoint <- function(a, b) {
  if (a < 0 && b > 0) {
    return(0)
  }

  mid <- (a + b) / 2
  width <- b - a
  mag <- floor(log10(width))

  for (m in mag + 5:-5) {
    step <- 10^m

    first <- floor(a / step) * step

    if (first + step > b) {
      next
    }

    # All multiples of step in [a, b]
    candidates <- seq(first, b, by = step)
    candidates <- candidates[candidates > a & candidates <= b]

    if (!is_null(candidates)) {
      # Pick the one closest to the midpoint
      return(candidates[which.min(abs(candidates - mid))])
    }
  }

  mid
}

diff1 <- function(x) {
  x[-1L] - x[-length(x)]
}

pseudo_inverse <- function(V, tol = 1e-9) {

  Vsvd <- svd(V)
  pos <- which(Vsvd$d > max(tol * Vsvd$d[1L], 0))

  if (is_null(pos)) {
    V[] <- 0
    return(V)
  }

  ss(Vsvd$v, j = pos) %*% ((1 / Vsvd$d[pos]) * t(ss(Vsvd$v, j = pos)))
}

# The branch path of a tree node, as a single label: the split conditions leading to it,
# joined with `&`. The conditions themselves can contain commas (a set of factor levels),
# so `&` rather than `,` keeps them readable. Returns NULL at the root, which has no
# branch path.
path_label <- function(path, collapse = " & ") {
  if (is_null(path)) {
    return(NULL)
  }

  if (identical(path[1L], "root")) {
    path <- path[-1L]
  }

  if (is_null(path)) {
    return(NULL)
  }

  paste(path, collapse = collapse)
}

# The ERF stored at a node of an <erf_tree> is a thin <erf>: enough to plot that node's
# dose-response and nothing else. It inherits from "erf", so functions that need the
# data, the weights, or the fitted nuisance models have to reject it explicitly.
check_full_erf <- function(x, arg = "erf") {
  if (inherits(x, "erf_node")) {
    arg::err("{.arg {arg}} is the ERF at a single node of an {.cls erf_tree}, which carries only what is needed to plot that node; supply the {.cls erf} returned by {.fn estimate_erf} instead")
  }

  invisible(x)
}

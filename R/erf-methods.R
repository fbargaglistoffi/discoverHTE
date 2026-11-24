#' Update and Display an Estimated ERF
#' @name update.erf
#'
#' @description
#' `update()` re-estimates an `<erf>`, and does so without refitting the nuisance models
#' when the bandwidth is the only thing that changes. `print()` displays the variables, the
#' sample size, the exposure range, and the nuisance method of an `<erf>`.
#'
#' @param object an `<erf>` object; the output of a call to [estimate_erf()].
#' @param x an `<erf>` object.
#' @param ... for `update()`, arguments to change in the original call to [estimate_erf()];
#'   supplying `bw` alone takes the fast path described in Details. For `print()`, ignored.
#'
#' @return
#' For `update()`, an `<erf>` object with the same components as the one supplied. For
#' `print()`, `x` invisibly.
#'
#' @details
#' When `bw` is the only argument supplied, the ERF is re-smoothed at the new bandwidth from
#' the nuisance estimates already stored in `object`. The nuisance models are not fit again,
#' which makes this far cheaper than a fresh call to [estimate_erf()] and is the usual
#' reason to call `update()` on an `<erf>`. The bandwidth is processed as it is by
#' [estimate_erf()], so a number, a function implementing a rule of thumb (such as
#' [bw.SJ()]), or `NULL` to select it by leave-one-out cross-validation are all allowed.
#'
#' Any other argument sends the call to [stats::update.default()], which re-evaluates the
#' original call with the new arguments in place. The nuisance models are fit again in that
#' case, so the cost is that of the original call.
#'
#' The ERF at a node of an `<erf_tree>`, as extracted by [`[[.erf_tree`], cannot be updated;
#' it carries the node's estimates but not the fitted nuisance models that re-smoothing
#' requires, and supplying one is an error.
#'
#' @seealso [estimate_erf()], [plot.erf()]
#'
#' @examples
#' set.seed(1234)
#' n <- 500
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' A <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)
#' Y <- 2 * A + X1 + X2 + rnorm(n)
#' d <- data.frame(Y, A, X1, X2)
#'
#' e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 25L)
#' e
#'
#' # Re-smooth at a wider bandwidth, reusing the fitted nuisance models
#' e_wide <- update(e, bw = 2 * e$bw)
#' e_wide$bw
#'
#' # The estimates change; the nuisance fit behind them does not
#' e$erf(1, vcov = FALSE)$est
#' e_wide$erf(1, vcov = FALSE)$est
#'
#' @exportS3Method print erf
print.erf <- function(x, ...) {
  cli::format_inline("An {.cls erf} object\n") |>
    cli::style_italic() |>
    cli::cat_line()

  nm <- get_formula_vars(x$formula)
  method <- x$method |>
    attr("method_name")

  cli::format_inline(" - Variables:\n",
                     "   + outcome: {nm$outcome}\n",
                     "   + exposure: {nm$exposure}\n",
                     "   + covariates: {nm$covariates}\n",
                     " - Sample size: {nrow(x$data)}\n",
                     " - Exposure range: {round(min(x$values), 4)} to {round(max(x$values), 4)}\n",
                     " - Nuisance method: {method}") |>
    cli::cat_line()

  invisible(x)
}

#' @rdname update.erf
#' @exportS3Method stats::update erf
update.erf <- function(object, ...) {
  # Re-smoothing needs the fitted nuisance models, which a node ERF does not carry.
  check_full_erf(object, "object")

  if (...length() == 1L && "bw" %in% ...names()) {
    call <- match.call()

    bw <- process_bw(...elt(1), A = object[["treat"]])

    # Only the bandwidth changes. Everything else the ERF was computed under has to be
    # supplied again, because estimate_erf_internal() defaults each of these to the
    # unclustered, zeroth-order, untrimmed case rather than reading them off `object`.
    # The clustering is rebuilt the way estimate_erf() builds it, from the clustering
    # stored on the analysis rows.
    cluster <- process_cluster(cluster = object[["cluster"]],
                               data = object[["cluster"]])

    # Pass to internal function
    erf <- estimate_erf_internal(fit = object[["fit"]],
                                 values = object[["values"]],
                                 bw = bw,
                                 cluster = cluster,
                                 verbose = FALSE,
                                 d = object[["d"]] %or% 0,
                                 trim = object[["trim"]] %or% 1)

    object[["erf"]] <- erf

    # The resolved bandwidth, not the argument: `bw` may be NULL (select by
    # cross-validation) or a function implementing a rule of thumb.
    object[["bw"]] <- attr(erf, "bw")

    # `[[<-` on a call deletes the element when the value is NULL, which errors when the
    # original call had no `bw` -- so `update(object, bw = NULL)` used to fail outright.
    object[["call"]] <- rlang::call_modify(object[["call"]], bw = call[["bw"]])

    return(object)
  }

  NextMethod("update")
}

process_nuissance <- function(method) {
  method <- force(method)

  arg::arg_is(method, "erf_nuisance")

  method
}

process_exposure <- function(formula, data) {
  f <- Formula::Formula(formula)

  l <- length(f)

  if (length(l) != 2L || l[1L] != 1L || l[2L] != 2L) {
    arg::err("{.arg formula} must be structured like {.code Y ~ A | X1 + X2}, where {.code Y} is the outcome, {.code A} is the exposure, and {.code X1} and {.code X2} are the covariates")
  }

  A <- model.frame(f, data = data, lhs = 0, rhs = 1)

  if (ncol(A) != 1L) {
    arg::err("the exposure must be a single variable")
  }

  A_name <- names(A)

  A <- A[[1L]]

  if (!is.numeric(A) || fnunique(A) <= 2L) {
    arg::err("the exposure must be a numeric variable with more than 2 unique values")
  }

  attr(A, "varname") <- A_name

  A
}

process_formula <- function(formula, data) {
  f <- Formula::Formula(formula)

  l <- length(f)

  if (length(l) != 2L || l[1L] != 1L || l[2L] != 2L) {
    arg::err("{.arg formula} must be structured like {.code Y ~ A | X1 + X2}, where {.code Y} is the outcome, {.code A} is the exposure, and {.code X1} and {.code X2} are the covariates")
  }

  # `na.action = na.pass` keeps rows with missing values so that `Y`, `A`, and `X`
  # stay full-length and row-aligned. With no missing data this is identical to the
  # default `na.omit` (no rows to drop); it matters for methods like
  # `method_negbinomcens()` where a missing outcome flags an interval-censored cell
  # that must be retained, not dropped.
  Y <- model.frame(f, data = data, lhs = 1, rhs = 0, na.action = na.pass)

  if (ncol(Y) != 1L) {
    arg::err("the outcome must be a single variable")
  }

  Y_name <- names(Y)

  A <- model.frame(f, data = data, lhs = 0, rhs = 1, na.action = na.pass)

  if (ncol(A) != 1L) {
    arg::err("the exposure must be a single variable")
  }

  A_name <- names(A)

  if (Y_name %in% A_name) {
    arg::err("the outcome cannot be the same as the exposure")
  }

  X <- model.frame(f, data = data, lhs = 0, rhs = 2, na.action = na.pass)

  if (Y_name %in% names(X)) {
    arg::err("the outcome cannot be in the covariates")
  }

  if (A_name %in% names(X)) {
    arg::err("the exposure cannot be in the covariates")
  }

  Y <- Y[[1L]]
  A <- A[[1L]]

  if (!is.numeric(A) || fnunique(A) <= 2L) {
    arg::err("the exposure must be a numeric variable with more than 2 unique values")
  }

  attr(A, "varname") <- A_name
  attr(Y, "varname") <- Y_name

  list(Y = Y,
       A = A,
       X = X)
}

get_formula_vars <- function(formula) {
  f <- Formula::Formula(formula)

  outcome <- formula(f, lhs = 1L, rhs = 0L) |>
    get_varnames()

  exposure <- formula(f, lhs = 0L, rhs = 1L) |>
    get_varnames()

  covariates <- formula(f, lhs = 0L, rhs = 2L) |>
    get_varnames()

  list(outcome = outcome,
       exposure = exposure,
       covariates = covariates)
}

process_range <- function(range = .95, n, A, w = NULL, strict = FALSE) {
  arg::arg_numeric(range)
  arg::arg_length(range, 1:2)


  arg::arg_count(n)
  arg::arg_gte(n, 2)

  if (length(range) > 1L) {
    range <- sort(range)

    range_t <- .range(A)

    if (strict) {
      if (any(range < range_t[1L]) || any(range > range_t[2L])) {
        r <- format(range_t, digits = 4L, drop0trailing = TRUE)
        arg::err("no values in {.arg range} can be outside the range of observed exposure values ({r[1L]} to {r[2L]})")
      }
    }
    else {
      range_diff <- range_t[2L] - range_t[1L]

      if (any(range < range_t[1L] - .1 * range_diff) ||
          any(range > range_t[2L] + .1 * range_diff)) {
        # arg::err("some values in {.arg range} are outside the range of observed exposure values")
        arg::wrn("some values in {.arg range} are outside the range of observed exposure values")
      }
    }
  }
  else if (range == 1) {
    range <- .range(A)
  }
  else if (all(range > 0 | range < 1)) {
    range <- .quantile(A, w = w, probs = c(1 - range, 1 + range) / 2)
  }
  else {
    arg::err("if {.arg range} is a single number, it must be a quantile (i.e., between 0 and 1)")
  }

  seq(range[1L], range[2L], length.out = n)
}

process_bw <- function(bw, A) {
  arg::when_not_null(bw,
                     arg::arg_or(
                       arg::arg_number,
                       arg::arg_function
                     ))

  if (is_null(bw)) {
    return(NULL)
  }

  if (is.function(bw)) {
    test_out <- try(bw(A), silent = TRUE)

    if (!is.numeric(test_out) || length(test_out) != 1L || anyNA(test_out)) {
      arg::err("the function supplied to {.arg bw} must return a single number")
    }
  }

  bw
}

process_crossfit <- function(crossfit) {
  arg::arg_or(crossfit,
              arg::arg_flag,
              arg::arg_count)

  if (isTRUE(crossfit)) {
    return(2L)
  }

  if (isFALSE(crossfit)) {
    return(1L)
  }

  as.integer(crossfit)
}

process_d <- function(d) {
  arg::arg_number(d)

  # Only the ERF and its first derivative are offered. The local polynomial order is tied
  # to `d` as p = d + 1, and each additional order needs its own bandwidth behavior to be
  # worth exposing.
  if (!d %in% c(0, 1)) {
    arg::err("{.arg d} must be {.val {0}} for the ERF itself or {.val {1}} for its first derivative")
  }

  as.integer(d)
}

process_trim <- function(trim) {
  arg::arg_number(trim)
  arg::arg_gt(trim, 0)
  arg::arg_lte(trim, 1)

  trim
}

process_weights <- function(weights, data) {
  if (is_null(weights)) {
    return(alloc(1.0, nrow(data)))
  }

  if (rlang::is_string(weights) && weights %in% names(data)) {
    weights <- data[[weights]]

    if (!is.numeric(weights)) {
      arg::err("the variable named in {.arg weights} must be numeric")
    }
  }

  if (!is.numeric(weights) || length(weights) != nrow(data)) {
    arg::err("{.arg weights} must be the name of a variable in {.arg data} or a numeric vector with a value for each unit in {.arg data}")
  }

  weights
}

process_cluster <- function(cluster, data) {
  if (is_null(cluster)) {
    return(NULL)
  }

  if (rlang::is_formula(cluster)) {
    cluster <- try({
      model.frame(cluster, data = data, na.action = na.pass)
    }, silent = TRUE)

    if (null_or_error(cluster)) {
      arg::err("clusters could not be extracted from the data; try supplying {.arg cluster} as a data frame or list")
    }
  }
  else {
    cluster <- qDF(cluster)
  }

  if (fnrow(cluster) != fnrow(data)) {
    arg::err("the number of observations in {.arg cluster} must equal that in the supplied dataset")
  }

  if (anyNA(cluster, recursive = TRUE)) {
    arg::err("{.val {NA}}s are not allowed in the clustering variable")
  }

  # A single clustering variable only. Clustering also determines the cross-fitting folds
  # (whole clusters are kept together), which requires one partition of the sample;
  # non-nested multiway clustering admits no such partition, so accepting it would give
  # cluster-robust standard errors that the fold structure silently contradicts.
  if (fncol(cluster) > 1L) {
    arg::err("only one clustering variable is allowed in {.arg cluster}, but {fncol(cluster)} were supplied")
  }

  # Drop unused levels: qF() keeps them, and the small-sample adjustment below counts
  # levels, so a cluster variable carrying unused levels (routine after subsetting a data
  # frame) would otherwise overstate the number of clusters and mis-scale every SE.
  cluster[] <- lapply(cluster, function(x) droplevels(qF(x)))

  #For small sample adjustment (setting cadjust = TRUE in vcovCL)
  g <- vapply(cluster, fnlevels, numeric(1L))

  list(cluster = cluster,
       adj = g / (g - 1))
}

# The clustering variable as a one-column data frame aligned to the rows of `data`,
# resolved before the nuisance models are fit so that it can shape the cross-fitting
# folds as well as the standard errors. `NULL` means unclustered.
resolve_cluster <- function(cluster, data) {
  if (is_null(cluster)) {
    return(NULL)
  }

  if (rlang::is_formula(cluster)) {
    cluster <- try({
      model.frame(cluster, data = data, na.action = na.pass)
    }, silent = TRUE)

    if (null_or_error(cluster)) {
      arg::err("clusters could not be extracted from the data; try supplying {.arg cluster} as a vector")
    }
  }
  else {
    cluster <- qDF(cluster)
  }

  if (fncol(cluster) > 1L) {
    arg::err("only one clustering variable is allowed in {.arg cluster}, but {fncol(cluster)} were supplied")
  }

  if (fnrow(cluster) != fnrow(data)) {
    arg::err("the number of observations in {.arg cluster} must equal that in the supplied dataset")
  }

  if (anyNA(cluster, recursive = TRUE)) {
    arg::err("{.val {NA}}s are not allowed in the clustering variable")
  }

  cluster
}

# Cross-fitting fold assignment. When a clustering variable is supplied, whole clusters
# are assigned to folds together, so no unit's nuisance prediction is trained on its own
# cluster-mates; otherwise folds are assigned over rows.
make_folds <- function(K, n, cluster = NULL) {
  if (is_null(cluster)) {
    return(sample(rep_len(seq_len(K), n)))
  }

  cluster <- as.character(cluster)
  lev <- unique(cluster)

  setNames(sample(rep_len(seq_len(K), length(lev))), lev)[cluster] |>
    unname()
}

process_modifiers <- function(modifiers = NULL, erf) {
  # `na.pass` keeps the modifier frame the same height as the fitted nuisance, which the
  # node subsets index into. Dropping incomplete rows here would leave a shorter logical
  # that R recycles against the fit, silently shifting which units each node holds. Units
  # missing a modifier are handled by the splitters, which send them to the right branch.
  if (is_null(modifiers)) {
    f <- Formula::Formula(erf$formula)
    return(model.frame(f, data = erf$data, lhs = 0, rhs = 2, na.action = na.pass))
  }

  f <- Formula::Formula(modifiers)

  l <- length(f)

  if (length(l) != 2L || l[1L] != 0L || l[2L] != 1L) {
    arg::err("{.arg modifiers} must be structured like `~ M1 + M2`, where `M1` and `M2` are the effect modifiers")
  }

  model.frame(f, data = erf$data, lhs = 0, rhs = 1, na.action = na.pass)
}

process_metric <- function(metric) {
  arg::match_arg(metric, c("L2", "L1", "AUC", "fANOVA"))
}

process_int_weight_method <- function(int_weight_method) {
  arg::match_arg(int_weight_method, c("none", "density", "gps"))
}

make_labels <- function(levels_in, levels_out, v) {
  if (length(levels_in) == 1L) {
    if (length(levels_out) == 1L) {
      return(c(sprintf("%s = %s", v, levels_in),
               sprintf("%s = %s", v, levels_out)))
    }
    else {
      return(c(sprintf("%s = %s", v, levels_in),
               sprintf("%s \u2260 %s", v, levels_in)))
    }
  }

  if (length(levels_out) == 1L) {
    return(c(sprintf("%s \u2260 %s", v, levels_out),
             sprintf("%s = %s", v, levels_out)))
  }

  c(sprintf("%s \u2208 {%s}", v, toString(levels_in)),
    sprintf("%s \u2208 {%s}", v, toString(levels_out)))
}

process_coords <- function(coords, cutoff, data, cluster) {
  if (is_null(coords)) {
    if (!is_null(cutoff)) {
      arg::err("{.arg cutoff} was supplied without {.arg coords}")
    }
    
    return(NULL)
  }
  
  if (is_null(cluster)) {
    arg::err("{.arg cluster} must be supplied when {.arg coords} is")
  }
  
  if (!is.character(coords) || length(coords) != 2L || !all(coords %in% names(data))) {
    arg::err("{.arg coords} must name the longitude and latitude variables in {.arg data}, in that order")
  }
  
  arg::arg_number(cutoff)
  arg::arg_gt(cutoff, 0)
  
  lon <- data[[coords[1L]]]
  lat <- data[[coords[2L]]]
  
  if (!is.numeric(lon) || !is.numeric(lat) || anyNA(lon) || anyNA(lat) || any(abs(lat) > 90)) {
    arg::err("the {.arg coords} variables must be non-missing decimal degrees, longitude first")
  }
  
  locs <- unique(data.frame(cl = as.character(cluster[[1L]]), lon = lon, lat = lat))
  
  if (anyDuplicated(locs$cl) > 0L) {
    arg::err("the longitude and latitude must not vary within a cluster")
  }
  
  cbind(lon = as.numeric(lon), lat = as.numeric(lat))
}


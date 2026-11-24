sim_count_data <- function(n_cty = 30, n_unit = 6, seed = 2) {
  set.seed(seed)

  ages <- c("65-74", "75-84", "85+")
  cty <- rep(sprintf("C%02d", seq_len(n_cty)), each = n_unit)
  unit <- paste0(cty, "_", ave(seq_along(cty), cty, FUN = seq_along))

  d <- do.call("rbind", lapply(seq_along(unit), function(i) {
    data.frame(unit = unit[i], cty = cty[i], age = ages,
               pop = round(runif(3, 2e3, 1e4)),
               pm = rnorm(1, 10, 2), z = rnorm(1))
  }))

  d$cnt <- rpois(nrow(d), d$pop * exp(-6 + 0.03 * d$pm + 0.1 * d$z))
  d
}

test_that("clustering supplied to estimate_erf reaches the unit level", {
  skip_if_not_installed("mgcv")

  d <- sim_count_data()
  m <- method_negbinomcens(age = "age", population = "pop", unit = "unit")

  e <- estimate_erf(cnt ~ pm | z, data = d, method = m,
                    cluster = ~cty, crossfit = 3, n = 15L)

  # The analysis rows are units, so the stored clustering must be too.
  expect_identical(length(e$treat), length(unique(d$unit)))
  expect_identical(nrow(e$cluster), length(e$treat))
  expect_identical(length(unique(as.character(e$cluster[[1L]]))),
                   length(unique(d$cty)))
})

test_that("cluster is no longer an argument of method_negbinomcens", {
  skip_if_not_installed("mgcv")

  d <- sim_count_data()
  m <- method_negbinomcens(age = "age", population = "pop", unit = "unit",
                           cluster = "cty")

  # `...` is funnelled into exposure/outcome, so this must be caught explicitly rather
  # than silently absorbed.
  expect_error(m(formula = cnt ~ pm | z, data = d, crossfit = 0L, n_grid = 20L),
               regexp = "no longer an argument", ignore.case = TRUE)
})

test_that("a cluster that varies within a unit is rejected", {
  skip_if_not_installed("mgcv")

  d <- sim_count_data()
  d$per_row <- factor(seq_len(nrow(d)))     # varies within every unit

  m <- method_negbinomcens(age = "age", population = "pop", unit = "unit")

  expect_error(estimate_erf(cnt ~ pm | z, data = d, method = m,
                            cluster = ~per_row, n = 15L),
               regexp = "must not vary within", ignore.case = TRUE)
})

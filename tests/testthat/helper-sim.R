# Small simulated dataset used across tests. `mod` is a genuine effect modifier of the
# exposure-response slope; `cl` is a clustering variable with `n / clsize` clusters.
sim_erf_data <- function(n = 400, clsize = 10, seed = 1234) {
  set.seed(seed)

  X1 <- rnorm(n)
  X2 <- rnorm(n)
  mod <- factor(rbinom(n, 1, 0.5))
  A <- 1 + 0.5 * X1 + rnorm(n)
  Y <- A * (1 + 1.5 * (mod == "1")) + X1 + X2 + rnorm(n)

  data.frame(Y, A, X1, X2, mod,
             cl = factor(rep(seq_len(ceiling(n / clsize)), each = clsize)[seq_len(n)]))
}

# Data whose ERF is known in closed form. `slope` makes it linear with that constant
# derivative; `quad = TRUE` makes it theta(a) = a^2, whose derivative is 2a and so varies
# across the exposure range. The covariates confound the exposure, so an unadjusted
# regression of Y on A does not recover the ERF.
sim_known_erf_data <- function(n = 2000, slope = 2, quad = FALSE, seed = 7) {
  set.seed(seed)

  X1 <- rnorm(n)
  X2 <- rnorm(n)
  A <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)

  Y <- {
    if (quad) A^2 + X1 + X2 + rnorm(n)
    else slope * A + X1 + X2 + rnorm(n)
  }

  data.frame(Y, A, X1, X2)
}

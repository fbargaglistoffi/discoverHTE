.simultaneous_conf_level <- function(vcov, conf_level = .95, maxiter = 1e5, ptol = 1e-4, ...) {
  if (conf_level == 0) {
    return(0)
  }

  arg::arg_supplied(vcov)
  arg::arg_cov(vcov)

  # Convert to correlation
  if (!all(check_if_zero(diag(vcov) - 1)) ||
      any(vcov < -1) || any(vcov > 1)) {
    zeros <- check_if_zero(diag(vcov))

    vcov <- vcov |>
      ss(!zeros, !zeros) |>
      cov2cor()
  }

  withCallingHandlers({
    z_crit <- mvtnorm::qmvnorm(conf_level,
                               tail = "both.tails",
                               corr = vcov,
                               maxiter = maxiter,
                               ptol = ptol,
                               ...)
  },
  error = function(e) {
    arg::err("there was an error computing simultaneous confidence intervals:\n{conditionMessage(e)}")
  })

  1 - 2 * pnorm(-z_crit$quantile)
}

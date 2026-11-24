#' Plot an Exposure-Response Function
#'
#' @description
#' Plots an ERF and its confidence band.
#'
#' @param x an `<erf>` object; the output of a call to [estimate_erf()].
#' @param conf_level the desired confidence level. Default is .95 for 95% confidence bands. Set to 0 to omit confidence bands.
#' @param simultaneous `logical`; whether the computed confidence bands should be simultaneous (`TRUE`) or pointwise (`FALSE`). Simultaneous (also known as uniform) bands cover the full line at the desired confidence level, whereas pointwise confidence bands only cover each point at the desired level. Default is `TRUE`.
#' @param density `logical`; whether to additionally plot the exposure density. Default is `TRUE`.
#' @param ... Ignored.
#'
#' @return
#' A `<ggplot>` object that can modified using functions in \pkg{ggplot2}.
#'
#' @details
#' `plot()` displays the ERF in a plot. The solid line corresponds to the ERF and the ribbon around it corresponds to its confidence band.
#'
#' @seealso [estimate_erf()]
#'
#' @examplesIf rlang::is_installed(c("ggside", "mvtnorm"))
#' set.seed(1234)
#' n <- 500
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' A <- 1 + 0.5 * X1 - 0.3 * X2 + rnorm(n)
#' Y <- 2 * A + X1 + X2 + rnorm(n)
#' d <- data.frame(Y, A, X1, X2)
#'
#' e <- estimate_erf(Y ~ A | X1 + X2, data = d, n = 25L)
#'
#' # Simultaneous bands and the exposure density, the defaults
#' plot(e)
#'
#' # Pointwise bands, and no density along the top
#' plot(e, simultaneous = FALSE, density = FALSE)
#'
#' # The ERF on its own
#' plot(e, conf_level = 0, density = FALSE)
#'
#' # The result is a <ggplot> and can be modified as one
#' plot(e, density = FALSE) +
#'   ggplot2::labs(title = "Estimated exposure-response function",
#'                 y = "Expected outcome")
#'
#' @exportS3Method plot erf
plot.erf <- function(x, conf_level = .95, simultaneous = TRUE, density = TRUE, ...) {
  arg::arg_number(conf_level)
  arg::arg_between(conf_level, c(0, 1), inclusive = c(TRUE, FALSE))
  inf <- conf_level > 0

  if (!inf) {
    simultaneous <- FALSE
  }

  arg::arg_flag(simultaneous)
  arg::arg_flag(density)

  est <- x$erf(x$values, vcov = inf)

  plot_dat <- data.frame(
    treat = x$values,
    est = est$est
  )

  p <- ggplot()

  if (inf) {
    if (simultaneous) {
      rlang::check_installed("mvtnorm")
      conf_level <- .simultaneous_conf_level(est$vcov, conf_level = conf_level)
    }

    z_crit <- abs(qnorm((1 - conf_level) / 2))

    plot_dat$se <- sqrt(diag(est$vcov))

    plot_dat$ci_lo <- plot_dat$est - z_crit * plot_dat$se
    plot_dat$ci_hi <- plot_dat$est + z_crit * plot_dat$se

    p <- p +
      geom_line(aes(y = .data$est, x = .data$treat),
                data = plot_dat) +
      geom_ribbon(aes(ymin = .data$ci_lo, ymax = .data$ci_hi,
                    x = .data$treat),
                  data = plot_dat,
                alpha = .2)
  }
  else {
    p <- p +
      geom_line(aes(y = .data$est, x = .data$treat),
                data = plot_dat)
  }

  p <- p +
    coord_cartesian(expand = c(left = TRUE, right = TRUE), default = TRUE) +
    theme_bw()  +
    labs(x = attr(x$treat, "varname"),
         y = "E[Y]")

  if (density) {
    rlang::check_installed("ggside")

    p <- p +
      ggside::ggside(x.pos = "top", collapse = "x") +
      ggside::geom_xsidedensity(data = data.frame(treat = x$treat),
                                aes(x = .data$treat),
                                fill = "white",
                                inherit.aes = FALSE) +
      ggside::theme_ggside_bw() +
      theme(ggside.axis.text.y = element_blank(),
            ggside.axis.ticks.y = element_blank(),
            ggside.panel.scale = .2)
  }

  p
}

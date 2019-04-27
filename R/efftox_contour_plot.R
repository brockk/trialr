#' @title Plot EffTox utility contours
#'
#' @description Plot EffTox utility contours. The probability of efficacy is
#' on the x-axis and toxicity on the y-axis. The zero-utility curve is plotted
#' bolder. The three "hinge points" are plotted as blue triangles. Optional
#' Prob(Efficacy) vs Prob(Toxicity) points can be added; these are shown as
#' red numerals, enumerated in the order provided.
#'
#' @param fit An instance of \code{\link{efftox_fit}}.
#' @param use_ggplot logical, TRUE to use ggplot2. Defaults to FALSE to use
#'  standard R graphics.
#' @param prob_eff vector of numbers between 0 and 1, containing the
#' efficacy probabilities of extra points to add to the plot as points,
#' e.g. the posterior mean efficacy probabilities of the doses under
#' investigation. Paired with prob_tox, thus they should be the same length.
#' Defaults to the values fitted by the model. Use NULL to supress.
#' @param prob_tox vector of numbers between 0 and 1, containing the
#' toxicity probabilities of extra points to add to the plot as points,
#' e.g. the posterior mean toxicity probabilities of the doses under
#' investigation. Paired with prob_eff, thus they should be the same length.
#' Defaults to the values fitted by the model. Use NULL to supress.
#' @param num_points integer for number of points to calculate on each curve.
#' The default is 1000 and this should be plenty.
#' @param util_vals A contour is plotted for each of these utility values.
#' The default is contours spaced by 0.2 between from -3 and 3,
#' i.e. \code{seq(-3, 3, by = 0.2)}.
#'
#' @return if \code{use_ggplot = TRUE}, an instance of \code{ggplot}; else no
#' object is returned. Omit assignment in either case to just view the plot.
#'
#' @export
#'
#' @examples
#' fit <- stan_efftox_demo(outcome_str = '1N 2E 3B')
#' efftox_contour_plot(fit)
#' title('EffTox utility contours')
#' # The same with ggplot2
#' efftox_contour_plot(fit, use_ggplot = TRUE) +
#'                     ggplot2::ggtitle('EffTox utility contours')
#'
#' @seealso
#' \code{\link{stan_efftox}}
efftox_contour_plot <- function(fit,
                                use_ggplot = FALSE,
                                prob_eff = fit$prob_eff,
                                prob_tox = fit$prob_tox,
                                num_points = 1000,
                                util_vals = seq(-3, 3, by = 0.2)) {
  eff_vals = seq(0, 1, length.out = num_points)
  dat <- fit$dat

  if(!is.null(prob_eff) & !is.null(prob_tox))
    if(length(prob_eff) != length(prob_tox))
      stop('prob_eff and prob_tox should be the same length')

  if(use_ggplot) {
    tox_vals = sapply(util_vals, function(u) efftox_get_tox(eff_vals, u, dat$p,
                                                            dat$eff0, dat$tox1))
    df = data.frame(eff_vals = rep(eff_vals, times = length(util_vals)),
                    tox_vals = as.numeric(tox_vals),
                    util_vals = rep(util_vals, each = length(eff_vals)))

    plt <- ggplot2::ggplot(df, ggplot2::aes(x = eff_vals, y = tox_vals,
                                            group = as.factor(util_vals))) +
      ggplot2::geom_line(size = 0.5, alpha = 0.25) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::xlab('Prob(Efficacy)') +
      ggplot2::ylab('Prob(Toxicity)')

    # Add neutral utility contour
    tox_vals = efftox_get_tox(eff_vals, 0, dat$p, dat$eff0, dat$tox1)
    df2 = data.frame(eff_vals, tox_vals, util_vals = 0)
    plt <- plt + ggplot2::geom_line(data = df2, size = 1)

    # Add hinge points
    df3 <- data.frame(prob_eff = c(dat$eff0, 1, dat$eff_star),
                      prob_tox = c(0, dat$tox1, dat$tox_star))
    plt <- plt + ggplot2::geom_point(data = df3, ggplot2::aes(x = prob_eff,
                                                              y = prob_tox,
                                                              group = 1),
                                     col = 'blue', shape = 24, size = 3)

    # Add provided eff & tox points
    if(!is.null(prob_eff) & !is.null(prob_tox)) {
      df4 <- data.frame(prob_eff, prob_tox, dl = 1:length(prob_eff))
      dl <- NULL
      plt <- plt + ggplot2::geom_text(data = df4, ggplot2::aes(x = prob_eff,
                                                               y = prob_tox,
                                                               group = 1,
                                                               label = dl),
                                      col = 'red', size = 4)
    }
    return(plt)
  } else {
    graphics::plot(NULL, ylim = c(0, 1), xlim = c(0, 1), ylab = 'Prob(Toxicity)',
                   xlab = 'Prob(Efficacy)')

    for(u in util_vals) {
      tox_vals = efftox_get_tox(eff_vals, u, dat$p, dat$eff0, dat$tox1)
      graphics::points(eff_vals, tox_vals, type = 'l', col = 'grey', lwd = 0.2)
    }

    # Add neutral utility contour
    tox_vals = efftox_get_tox(eff_vals, 0, dat$p, dat$eff0, dat$tox1)
    graphics::points(eff_vals, tox_vals, type = 'l', col = 'black', lwd = 2)

    # Add hinge points
    graphics::points(dat$eff0, 0, col = 'blue', pch = 2)
    graphics::points(1, dat$tox1, col = 'blue', pch = 2)
    graphics::points(dat$eff_star, dat$tox_star, col = 'blue', pch = 2)

    # # Add provided eff & tox points
    if(!is.null(prob_eff) & !is.null(prob_tox)) {
      graphics::points(prob_eff, prob_tox, col = 'red',
                       pch = as.character(1:length(prob_eff)))
    }
  }
}

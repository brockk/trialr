
#' @title Plot densities of EffTox dose utilities
#'
#' @description Plot densities of EffTox dose utilities. Optionally plot only a
#' subset of the doses by specifying the \code{doses} parameter. This function
#' requires ggplot2 be installed.
#'
#' @param fit An instance of \code{efftox_fit}.
#' @param doses optional, vector of integer dose-levels to plot. E.g. to plot
#' only dose-levels 1, 2 & 3 (and suppress the plotting of any other doses), use
#' \code{doses = 1:3}
#'
#' @return an instance of \code{ggplot}. Omit assignment to just view the plot.
#'
#' @export
#'
#' @note This function requires that ggplot2 be installed.
#'
#' @examples
#' fit <- stan_efftox_demo('1N 2E 3B')
#' efftox_utility_density_plot(fit) + ggplot2::ggtitle('My doses')  # Too busy?
#' # Specify subset of doses to make plot less cluttered
#' efftox_utility_density_plot(fit, doses = 1:3) + ggplot2::ggtitle('My doses')
#'
efftox_utility_density_plot <- function(fit, doses = NULL) {
  if(!('ggplot2' %in% utils::installed.packages()))
    stop('This function requires ggplot2 be installed.')

  u <- rstan::extract(fit$fit, par = 'utility')[[1]]
  df <- data.frame(Utility = as.numeric(u),
                   D = rep(fit$dose_indices, each = nrow(u))
  )
  df$Dose = factor(df$D)
  if(!is.null(doses))
    df = df[df$D %in% doses, ]
  Dose <- Utility <- NULL
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Utility, group = Dose,
                                        colour = Dose)) +
    ggplot2::geom_density()
  return(p)
}

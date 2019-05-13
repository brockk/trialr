

#' Class to hold the elements of a single dose-finding analysis residing in a
#' pathway of analyses.
#'
#' A pathway in a dose-finding trial is a series of successive analyses. For
#' instance, the model will likely be fit to all of the outcomes observed
#' at the end of the first cohort, the second cohort, etc. This class holds the
#' elements reflecting the analysis, and the place of this analysis in the
#' pathway.
#'
#' @name dose_finding_path_node-class
#' @aliases dose_finding_path_node
#' @docType class
#'
#' @param node_id An integer representing the id of this node in a pathway.
#' @param parent_node_id An integer representing the id of this node's parent in
#' the pathway.
#' @param depth An integer representing the depth of this node in the pathway,
#' where the root has depth 0.
#' @param outcomes A string representing the outcomes observed at the time of
#' analysis. See \code{\link{df_parse_outcomes}} for a description of syntax and
#' examples.
#' @param next_dose An integer representing the dose recommended by the model
#' for the next patient or cohort of patients.
#' @param fit Object obtained from fitting the dose-finding model to
#' \code{outcomes}.
#' @param parent_fit Object obtained from fitting the dose-finding model to
#' the outcomes of the parent node. Comparing to \code{fit} will oten be
#' valuable.
#'
#' @return Instance of class \code{dose_finding_path_node}
#' @export
#'
#' @examples
#' \dontrun{
#' parent_outcomes <- '1NNN'
#' outcomes <- '1NNN 2NNT'
#' target <- 0.25
#' skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
#' parent_fit <- stan_crm(outcome_str = parent_outcomes, skeleton = skeleton,
#'                        target = target, model = 'empiric', beta_sd = 1)
#' fit <- stan_crm(outcome_str = outcomes, skeleton = skeleton,
#'                 target = target, model = 'empiric', beta_sd = 1)
#' dose_finding_path_node(node_id = 2, parent_node_id = 1, depth = 1,
#'                        outcomes = outcomes, next_dose = fit$recommended_dose,
#'                        fit = fit, parent_fit = parent_fit)
#' }
dose_finding_path_node <- function(node_id, parent_node_id, depth, outcomes,
                                   next_dose, fit, parent_fit) {
  x <- list(.node = node_id,
            .parent = parent_node_id,
            .depth = depth,
            outcomes = outcomes,
            next_dose = next_dose,
            fit = fit,
            parent_fit = parent_fit)
  class(x) <- c("dose_finding_path_node", "list")
  x
}


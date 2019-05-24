#' Spread the information in dose_finding_paths object to a wide data.frame format.
#'
#' @param df Optional \code{data.frame} like that returned by
#' as_tibble(dose_finding_paths). Columns .depth, .node, .parent are required.
#' All other columns are spread with a suffix reflecting depth.
#' @param dose_finding_paths Optional instance of dose_finding_paths. Required
#' if `df` is null.
#' @param max_depth integer, maximum depth of paths to traverse.
#'
#' @return A data.frame
#'
#' @importFrom tibble as_tibble
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter select rename full_join rename_at vars starts_with
#'
#' @export
#'
#' @examples
#' \dontrun{
#' target <- 0.25
#' skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
#' paths <- crm_dtps(skeleton = skeleton, target = target, model = 'empiric',
#'                   cohort_sizes = c(1, 1), next_dose = 3, beta_sd = 1)
#' spread_paths(dose_finding_paths = paths)
#'
#' df <- as_tibble(paths)
#' spread_paths(df)
#' spread_paths(df %>% select(-fit, -parent_fit, -dose_index))
#' }
spread_paths <- function(df = NULL,
                         dose_finding_paths = NULL,
                         max_depth = NULL) {

  if(is.null(df) & is.null(dose_finding_paths))
    stop('Specify either df or dose_finding_paths.')
  if(is.null(df) & !is.null(dose_finding_paths))
    df <- as_tibble(dose_finding_paths)
  if(is.null(max_depth)) max_depth <- max(df$.depth)
  if(!all(c('.depth', '.node', '.parent') %in% colnames(df)))
    stop("Columns '.depth', '.node' and '.parent' are required.")

  .depth <- .parent <- .node <- Node <- NULL

  depth = 0
  wide_df <- df %>%
    filter(.depth == depth) %>%
    select(-.parent, -.depth) %>%
    rename_at(vars(-starts_with(".")), function(x) paste0(x, depth)) %>%
    rename(Node = .node)

  for(depth in 1:max_depth) {
    sub_df <- df %>%
      filter(.depth == depth) %>%
      rename_at(vars(-starts_with(".")), function(x) paste0(x, depth)) %>%
      select(-.depth)

    wide_df <- wide_df %>%
      full_join(sub_df,
                by = c('Node' = '.parent'),
                suffix = paste0(".", c(depth - 1, depth))) %>%
      select(-Node) %>%
      rename(Node = .node)
  }

  wide_df %>% select(-Node)
}



# The dose_finding_paths class itself is simply a list of dose_finding_path_node
# objects. It is so simple, it does not need documentation.
# The dose_finding_path_node object is documented.


#' Cast \code{dose_finding_paths} object to \code{\link[tibble]{tibble}}.
#'
#' @param x Object of class \code{dose_finding_paths}.
#' @param ... Extra args passed onwards.
#'
#' @return Object of class \code{\link[tibble]{tibble}}
#'
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom purrr map_chr
#' @export
as_tibble.dose_finding_paths <- function(x, ...) {
  fit <- NULL
  tibble(
    .node = map_dbl(x, '.node'),
    .parent = map_dbl(x, '.parent'),
    .depth = map_dbl(x, '.depth'),
    outcomes = map_chr(x, 'outcomes'),
    next_dose = map_dbl(x, 'next_dose'),
    fit = map(x, 'fit'),
    parent_fit = map(x, 'parent_fit'),
    dose_index = map(fit, 'dose_indices'),
    ...
  )
}

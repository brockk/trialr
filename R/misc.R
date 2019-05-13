
#' Get index of element in vector with value closest to a target
#'
#' @param vector Identify element in this numeric vector
#' @param target numeric target
#'
#' @return an integer indexing \code{vector}
#' @export
#'
#' @examples
#' closest_to_target(c(0.1, 0.2, 0.3), 0.05)  # 1
#' closest_to_target(c(0.1, 0.2, 0.3), 0.22)  # 2
#' closest_to_target(c(0.1, 0.2, 0.3), -0.05) # 1
#' closest_to_target(c(0.1, 0.2, 0.3), 8) # 3
closest_to_target <- function(vector, target) {
  which.min(abs(vector - target))
}

.entropy <- function(probs) {
  gt0 <- probs > 0
  -sum(probs[gt0] * log(probs[gt0]))
}

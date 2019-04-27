

.entropy <- function(probs) {
  gt0 <- probs > 0
  -sum(probs[gt0] * log(probs[gt0]))
}

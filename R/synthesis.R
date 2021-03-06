#' Synthesis operator.
#'
#' Compute the synthesis operator for coefficient coeff.
#'
#' @export synthesis
#' @param coeff Transform coefficients.
#' @param tf Frame coefficients.
#' @return \code{y} Synthesis signal.
#' @seealso \code{\link{analysis}}, \code{\link{tight_frame}}

synthesis <- function(coeff, tf) {
  y <- t(coeff) %*% tf
  return(as.vector(y))
}


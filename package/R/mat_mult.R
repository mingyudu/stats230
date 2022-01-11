#' Multiply two square matrices and a vector
#'
#' @param a a square matrix
#' @param b another square matrix
#' @param x a vector
#' @param opt argument indicating the order of multiplication
#' @return The multiplication of a, b, and x
#' @example
#' A <- matrix(2, 3, 3)
#' B <- matrix(4, 3, 3)
#' X <- matrix(5, 3, 1)
#' mat_mult(A, B, X, 1)
#' mat_mult(A, B, X, 2)
mat_mult <- function(a, b, x, opt){
  if(opt == 1){
    m <- a %*% b
    return(m %*% x)
  }
  if(opt == 2){
    m <- b %*% x
    return(a %*% m)
  }
}

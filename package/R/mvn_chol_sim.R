#' Simulate from a multivariate normal distribution using Cholesky decomposition
#'
#' @param n the number of realizations.
#' @param mu a vector giving the means of the variables.
#' @param sigma a positive definite symmetric matrix specifying the covariance matrix of the variables.
#' @return n realizations from the specified multivariate normal distribution.
#'
#' @import MASS
#' @export
#'
#' @examples
#' mu <- c(1,2)
#' sigma <- matrix(c(10,3,3,2),2,2)
#' mvn_chol_sim(10,mu,sigma)
mvn_chol_sim <- function(n,mu,sigma){
  l <- chol(sigma)
  z <- mvrnorm(n,rep(0,length(mu)),diag(1,length(mu)))
  x <- mu + l%*%t(z)
  return(x)
}

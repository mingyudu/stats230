#' OLS estimation of regression coefficients via singular value decomposition
#'
#' @param x design matrix of predictor variables
#' @param y vector of response variable
#' @return Estimated coefficients of linear regression model
ols_svd <- function(x, y){
  svd <- svd(x)
  v <- svd$v
  u <- svd$u
  sigma_inv <- 1/svd$d
  vec1 <- t(u) %*% y
  vec2 <- sigma_inv * vec1
  vec3 <- v %*% vec2
  return(vec3)
}

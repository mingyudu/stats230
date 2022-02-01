#' OLS estimation of regression coefficients via QR decomposition
#'
#' @param x design matrix of predictor varialbes
#' @param y vector of response variable
#' @return Estimated coefficients of linear regression model
ols_qr <- function(x, y){
  qr <- qr(x)
  q <- qr.Q(qr)
  r <- qr.R(qr)
  vec <- t(q) %*% y
  beta <- solve(r) %*% vec
  return(beta)
}

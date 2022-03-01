#' Logistic regression using gradient descent or Newton-Raphson method
#'
#' @param x design matrix of predictor variables(including the intercept term).
#' @param y vector of binary response variable.
#' @param threshold threshold value to test for convergence.
#' @param max_iter maximum number of iterations.
#' @param method indicator of method; if method=1, use gradient descent; if method=2, use Newton-Raphson method.
#' @param alpha learning rate for gradient descent.
#' @return MLEs of regression coefficients, their corresponding asymptotic confidence intervals, and a vector of the log-likelihoods “visited” by the chosen optimization option.
logistic.reg <- function(x, y, threshold, max_iter, method, alpha = 0.1){
  n <- length(y)
  beta <- rep(0, dim(x)[2])
  iter <- 0
  diff <- 10000
  llk <- NULL
  if(method==1){ # Gradient descent
    while (diff > threshold) {
      p <- as.vector(1/(1+exp(-x %*% beta)))

      # log-likelihood
      llk <- c(llk, sum(y*log(p) + (1-y)*log(1-p)))

      # update beta
      beta_diff <- (alpha/length(y))*(t(x) %*% (y-p))
      beta <- beta + beta_diff

      # calculate diff
      diff <- sum(beta_diff^2)
      iter <- iter + 1
      if(iter > max_iter){
        stop("The algorithm is not converging.")
      }
    }
    # confidence interval
    p <- as.vector(1/(1+exp(-x %*% beta)))
    w <- diag(p*(1-p))
    se <- sqrt(diag(solve(t(x) %*% w %*% x)))
    upper <- beta + 1.96 * se
    lower <- beta - 1.96 * se

    return(list(beta = round(beta, 3), llk = llk, se = se,
                conf.upper = upper, conf.lower = lower))
  }
  else{ # Newton-Raphson
    while (diff > threshold) {
      p <- as.vector(1/(1+exp(-x %*% beta)))

      # log-likelihood
      llk <- c(llk, sum(y*log(p) + (1-y)*log(1-p)))

      # update beta
      grad <- t(x) %*% (y-p)
      w <- diag(p*(1-p))
      hessian <- t(x) %*% w %*% x
      beta_diff <- solve(hessian) %*% grad
      beta <- beta + beta_diff

      # calculate diff
      diff <- sum(beta_diff^2)
      iter <- iter + 1
      if(iter > max_iter){
        stop("The algorithm is not converging.")
      }
    }
    # confidence interval
    p <- as.vector(1/(1+exp(-x %*% beta)))
    w <- diag(p*(1-p))
    se <- sqrt(diag(solve(t(x) %*% w %*% x)))
    upper <- beta + 1.96 * se
    lower <- beta - 1.96 * se

    return(list(beta = round(beta, 3), llk = llk, se = se,
                conf.upper = upper, conf.lower = lower))
  }
}

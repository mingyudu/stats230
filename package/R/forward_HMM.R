#' Compute the forward probabilities in hidden Markov model
#'
#' @param iniDist initial distribution of the hidden Markov model.
#' @param transProb transition probabilities of the hidden Markov model.
#' @param emisProb emission probabilities of the hidden Markov model.
#' @param y.lst a sequence of observed states.
#' @return a matrix of the forward probabilities. The first dimension is the state and the second is time.
forward_HMM <- function(iniDist, transProb, emisProb, y.lst){
  iniDist[is.na(iniDist)] <- 0
  transProb[is.na(transProb)] <- 0
  emisProb[is.na(emisProb)] <- 0
  N <- length(y.lst)
  S <- dim(transProb)[1]
  a <- array(NA, c(S, N))
  # Initialize a_1(i) = v(i)e(y_1|i) for i=1,...,S
  for (s in 1:S) {
    a[s, 1] <- log(iniDist[s] * emisProb[s, y.lst[1]])
  }
  # Iteration
  for (t in 2:N) {
    for (i in 1:S) {
      logsum = - Inf
      for (j in 1:S) {
        temp <- a[j, t-1] + log(transProb[j, i])
        # dynamic rescaling
        if (temp > - Inf){
          logsum <- temp + log(1 + exp(logsum - temp))
        }
      }
      a[i, t] <- log(emisProb[i, y.lst[t]]) + logsum
    }
  }
  return(a)
}

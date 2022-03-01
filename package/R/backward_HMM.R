#' Compute the backward probabilities in hidden Markov model
#'
#' @param transProb transition probabilities of the hidden Markov model.
#' @param emisProb emission probabilities of the hidden Markov model.
#' @param y.lst a sequence of observed states.
#' @return a matrix of the backward probabilities. The first dimension is the state and the second is time.
backward_HMM <- function(transProb, emisProb, y.lst){
  transProb[is.na(transProb)] <- 0
  emisProb[is.na(emisProb)] <- 0
  N <- length(y.lst)
  S <- dim(transProb)[1]
  b <- array(NA, c(S, N))
  # Calculate log-likelihood first, then do exponential transformation
  # Initialize b_n(i) = 1 for i=1,..,S
  for (s in 1:S) {
    b[s, N] <- log(1)
  }
  # Iteration
  for (t in (N-1):1) {
    for (i in 1:S) {
      logsum <- - Inf
      for (j in 1:S) {
        temp <- b[j, t+1] +
          log(transProb[i, j] * emisProb[j, y.lst[t+1]])
        # dynamic rescaling
        if (temp > - Inf){
          logsum = temp + log(1 + exp(logsum - temp))
        }
      }
      b[i, t] <- logsum
    }
  }
  return(b)
}

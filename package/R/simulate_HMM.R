#' Simulate hidden and observed states for a hidden Markov model
#'
#' @param iniDist initial distribution of the hidden Markov model
#' @param transProb transition probabilities of the hidden Markov model.
#' @param emisProb emission probabilities of the hidden Markov model.
#' @param hidden vector of hidden states.
#' @param observed vector of observed states.
#' @param length length of the simulated sequence of the hidden and observed states.
#' @param seed seed number to generate random number.
#' @return a sequence of hidden and observed states.
simulate_HMM <- function(iniDist, transProb, emisProb,
                         hidden, observed, length, seed){
  set.seed(seed)
  x.lst <- NULL
  y.lst <- NULL
  x.lst <- c(x.lst, sample(hidden, 1, prob = iniDist))
  for (i in 2:length) {
    x <- sample(hidden, 1, prob = transProb[x.lst[i-1],])
    x.lst <- c(x.lst, x)
  }
  for (i in 1:length) {
    y <- sample(observed, 1, prob = emisProb[x.lst[i],])
    y.lst <- c(y.lst, y)
  }
  return(list(x.lst=x.lst, y.lst=y.lst))
}

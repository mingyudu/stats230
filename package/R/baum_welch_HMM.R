#' Update the parameters of the hidden Markov model
#'
#' @param iniDist initial distribution of the hidden Markov model.
#' @param transProb transition probabilities of the hidden Markov model.
#' @param emisProb emission probabilities of the hidden Markov model.
#' @param y.lst a sequence of observed states.
#' @return updated initial distribution, transition probabilities, and emission probabilities of HMM
update_param_HMM <- function(iniDist, transProb, emisProb, y.lst){
  N <- length(y.lst)
  S <- dim(transProb)[1]
  Ob <- dim(emisProb)[2]
  transProb_new <- array(0, dim = dim(transProb))
  emisProb_new <- array(0, dim = dim(emisProb))
  emisProb_new[1,] <- rep(1/6, 6)
  iniDist_new <- rep(0, length(iniDist))

  a <- forward_HMM(iniDist, transProb, emisProb, y.lst)
  b <- backward_HMM(transProb, emisProb, y.lst)

  # calculate log(P(Y))
  py <- a[1, N]
  py2 <- a[2, N]
  if (py2 > - Inf){
    py <- py2 + log(1 + exp(py - py2))
  }

  # update transition matrix
  for (i in 1:S) {
    for (j in 1:S) {
      temp <- a[i, 1] + log(transProb[i, j]) +
        log(emisProb[j, y.lst[2]]) + b[j, 2]
      for (t in 2:(N-1)) {
        # calculate g_t(i,j)
        tempt <- a[i, t] + log(transProb[i, j]) +
          log(emisProb[j, y.lst[t+1]]) + b[j, t+1]
        if (tempt > - Inf){
          temp <- tempt + log(1 + exp(temp - tempt))
        }
      }
      temp <- exp(temp - py)
      transProb_new[i, j] <- temp
    }
  }

  # update emission matrix
  # the first row is fixed to 1/6 so i starts from 2
  for (i in 2:S) {
    for (m in 1:Ob) {
      temp <- - Inf
      for (t in 1:N) {
        if (m == y.lst[t]){
          tempe <- a[i, t] + b[i, t]
          if (tempe > - Inf){
            temp <- tempe + log(1 + exp(temp - tempe))
          }
        }
      }
      temp <- exp(temp - py)
      emisProb_new[i, m] <- temp
    }
  }

  # update initial distribution
  for (i in 1:S){
    temp <- - Inf
    tempi <- a[i, 1] + b[i, 1]
    if(tempi > - Inf){
      temp <- tempi + log(1 + exp(temp - tempi))
    }
    temp <- exp(temp - py)
    iniDist_new[i] <- temp
  }

  return(list(iniDist = iniDist_new, transProb = transProb_new,
              emisProb = emisProb_new))
}

#' Estimate the parameters in hidden Markov model using Baum-Welch algorithm
#'
#' @param iniDist initial distribution of the hidden Markov model.
#' @param transProb transition probabilities of the hidden Markov model.
#' @param emisProb emission probabilities of the hidden Markov model.
#' @param y.lst a sequence of observed states.
#' @param max_iter maximum number of iteration
#' @param delta threshold to test for convergence
#' @return updated initial distribution, transition probabilities, and emission probabilities of HMM
baum_welch_HMM <- function(iniDist, transProb, emisProb, y.lst,
                           max_iter=100, delta=1e-10){
  iniDist_temp <- iniDist
  iniDist_temp[is.na(iniDist_temp)] <- 0
  transProb_temp <- transProb
  transProb_temp[is.na(transProb_temp)] <- 0
  emisProb_temp <- emisProb
  emisProb_temp[is.na(emisProb_temp)] <- 0

  diff <- NULL
  for (i in 1:max_iter) {
    updatedHMM <- update_param_HMM(iniDist_temp, transProb_temp,
                                   emisProb_temp, y.lst)
    I <- updatedHMM$iniDist
    T <- updatedHMM$transProb
    E <- updatedHMM$emisProb

    I <- I/sum(I)
    T <- T/apply(T,1,sum)
    E <- E/apply(E,1,sum)

    d <- sqrt(sum((transProb_temp - T)^2)) + sqrt(sum((emisProb_temp[2,] - E[2,])^2)) +
      sqrt(sum((iniDist_temp - I)^2))

    diff <- c(diff, d)
    transProb_temp <- T
    emisProb_temp <- E
    iniDist_temp <- I
    if(d < delta){
      break
    }
  }

  iniDist_temp[is.na(iniDist_temp)] <- NA
  transProb_temp[is.na(transProb_temp)] <- NA
  emisProb_temp[is.na(emisProb_temp)] <- NA

  return(list(iniDist = iniDist_temp, transProb = transProb_temp,
              emisProb = emisProb_temp, diff = diff))
}

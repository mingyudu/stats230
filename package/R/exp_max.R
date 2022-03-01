#' Expectation-Maximization algorithm for the ABO blood type problem
#'
#' @param nA observed number of individuals with phenotype A.
#' @param nAB observed number of individuals with phenotype AB.
#' @param nB observed number of individuals with phenotype B.
#' @param nO observed number of individuals with phenotype O.
#' @param pA initial value of allele frequency probability of A.
#' @param pB initial value of allele frequency probability of B.
#' @return allele frequency estimates of A, B, O
exp_max <- function(nA, nAB, nB, nO, pA, pB){
  pO <- 1 - pA - pB
  iter <- 0
  diff <- 1
  while (sum(diff > 1e-5) > 0) {
    # E-step
    nAA <- nA * pA / (pA + 2 * pO)
    nAO <- nA * 2 * pO / (pA + 2 * pO)
    nBB <- nB * pB / (pB + 2 * pO)
    nBO <- nB * 2 * pO / (pB + 2 * pO)

    # M-step
    pA_new <- (2 * nAA + nAO + nAB) / (2 * (nA + nB + nO + nAB))
    pB_new <- (2 * nBB + nBO + nAB) / (2 * (nA + nB + nO + nAB))
    pO_new <- 1 - pA_new - pB_new

    # calculate diff
    diff <- abs(c(pA, pB, pO) - c(pA_new, pB_new, pO_new))
    iter <- iter + 1
    pA <- pA_new
    pB <- pB_new
    pO <- pO_new
  }
  return(list(pA = pA, pB = pB, pO = pO))
}

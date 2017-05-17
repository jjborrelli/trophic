#' Likelihood Estimation for Cascade Food Web Model
#'
#' @param A The food web adjacency matrix.
#'
#' @return The likelihood that the given adjacency matrix was produced by the Cascade model.
#' @export
#'
#' @examples
cascade_likelihood<- function(A){

  L <- sum(A)
  S <- nrow(A)
  CONN <- L/(S*(S-1))
  p <- (2*CONN*S)/(S-1)

  if(sum(A[lower.tri(A)]) == 0){

    P.A_p <- (p^L)*(1-p)^((S*(S-1))/2-L)
    return(log10(P.A_p))

  }else{

    N <- A
    N[lower.tri(N)] <- 0
    L <- sum(N)
    S <- nrow(N)
    CONN <- L/(S*(S-1))
    p <- (2*CONN*S)/(S-1)
    P.A_p <- (p^L)*(1-p)^((S*(S-1))/2-L)

    K <- A
    K[upper.tri(K)] <- 0
    L <- sum(K)
    S <- nrow(K)
    q <- L/(S*(S-1))
    P.K_q <- q^L*(1-q)^(S^2-L)

    L.pq_A <- P.A_p * P.K_q

  }

  return(log10(L.pq_A))
}

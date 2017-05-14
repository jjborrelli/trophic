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
  }
}

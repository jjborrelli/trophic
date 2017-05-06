#' Anarchy
#'
#' @param S Number of species in the community.
#' @param C Average path length used to get connectance.
#'
#' @return An adjacency matrix of a random network where any two species are connected with probability \code{C/S}.
#' @export
#'
#' @examples
#' anarchy(20, 1.86)
anarchy <- function(S, C){
  p <- C/S
  m <- matrix(rbinom(S^2, 1, p), nrow = S, ncol = S)
  return(m)
}


#' Democracy
#'
#' @param S Number of species in the community.
#' @param C Average path length used to get connectance.
#'
#' @return An adjacency matrix of a random network where any two species are connected with probability \code{C/S}. this model differs from \code{anarchy} by not allowing any trophic loops.
#' @export
#'
#' @examples
#' democracy(20, 1.86)
democracy <- function(S, C){
  p <- C/S
  a <- matrix(0, nrow = S, ncol = S)
  a[upper.tri(a)] <- rbinom(sum(upper.tri(a)), 1, p)
  return(a)
}


#' Cascade Food Web Model
#'
#' @param S Number of species in the community.
#' @param C Average path length used to get connectance.
#'
#' @return The adjacency matrix of a cascade model food web. Species are arranged along a niche axis and may consume any species with a niche value less than their own with a fixed probability of \code{2CS/(S-1)}, creating a triangular matrix.
#' @export
#'
#' @examples
#' cascade(20, .1)
cascade <- function(S, C){
  a <- matrix(0, nrow = S, ncol = S)
  p <- (2*C*S)/(S-1)
  for(i in 1:(S-1)){
    a[i, (i+1):S] <- rbinom(length((i+1):S), 1, p)
  }
  return(a)
}

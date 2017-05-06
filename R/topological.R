#' Anarchy
#'
#' @param S Number of species in the community.
#' @param C Average path length used to get connectance.
#'
#' @return An adjacency matrix of a random network where any two species are connected with probability \code{C/S}.
#' @export
#'
#' @section Reference:
#'Cohen, J. E., and C. M. Newman. 1985. A stochastic theory of community food webs: I. Models and aggregated data. Proceedings of the Royal Society B 224:421–448.
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
#' @section Reference:
#'Cohen, J. E., and C. M. Newman. 1985. A stochastic theory of community food webs: I. Models and aggregated data. Proceedings of the Royal Society B 224:421–448.
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
#' @param C The connectance, or fraction of realized links in the food web.
#'
#' @return The adjacency matrix of a cascade model food web. Species are arranged along a niche axis and may consume any species with a niche value less than their own with a fixed probability of \code{2CS/(S-1)}, creating a triangular matrix.
#' @export
#'
#' @section Reference:
#'Cohen, J. E., and C. M. Newman. 1985. A stochastic theory of community food webs: I. Models and aggregated data. Proceedings of the Royal Society B 224:421–448.
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


#' Cascade Food Web Model Type II
#'
#' @param S Number of species in the community.
#' @param C Average path length used to get connectance.
#'
#' @return The adjacency matrix of a cascade model food web where species higher on the niche axis feed on those lower with a fixed probability \code{C/S}.
#' @export
#'
#' @section Reference:
#'Cohen, J. E., and C. M. Newman. 1985. A stochastic theory of community food webs: I. Models and aggregated data. Proceedings of the Royal Society B 224:421–448.
#'
#' @examples
#' cascade_II(20, 1.86)
cascade_II <- function(S, C){
  a <- matrix(0, nrow = S, ncol = S)
  p <- C/S
  for(i in 1:(S-1)){
    a[i, (i+1):S] <- rbinom(length((i+1):S), 1, p)
  }
  return(a)
}


#' Generalized Cascade Food Web Model
#'
#' @param S Number of species in the community.
#' @param C The connectance, or fraction of realized links in the food web.
#'
#' @return The adjacency matrix of a cascade model food web where species higher on the niche axis feed on those lower with a fixed probability \code{1/(2 * C) - 1}.
#' @export
#'
#' @section Reference:
#' Stouffer, D. B., J. Camacho, R. Guimera, C. A. Ng, and L. A. N. Amaral. 2005. Quantitative patterns in the structure of model and empirical food webs. Ecology 86:1301–1311.
#'
#' @examples
#' gen_cascade(20, .1)
gen_cascade <- function(S, C){
  n.i <- sort(runif(S))
  B <- 1/(2 * C) - 1
  a <- matrix(0, nrow = S, ncol = S)
  for(i in 1:length(n.i)){
    a[n.i < n.i[i], i] <- rbinom(sum(n.i < n.i[i]), 1, rbeta(sum(n.i < n.i[i]),1,((1/(2*C))-1)))
  }
  return(a)
}

################################################
################################################
##### Visualize adjacency

#' Adjacency Matrix Visualization
#'
#' @param a Adjacency matrix
#'
#' @return A plot of the adjacency matrix where black cells indicate interactions
#' @export
#'
#' @examples
#' n <- niche(30, .2)
#' avis(n)
avis <- function(a){
  bipartite::visweb(a, type = "none", clear = F)
}

################################################
################################################
##### Visualize network

#' Adjacency Matrix Visualization
#'
#' @param a Adjacency matrix
#'
#' @return A plot of the adjacency matrix where black cells indicate interactions
#' @export
#'
#' @examples
#' n <- niche(30, .2)
#' avis(n)
netvis_trophic <- function(a){
  net <- igraph::graph.adjacency(a)
  tl <- NetIndices::TrophInd(a)
  lay <- layout.kamada.kawai(net)
  lay[,2] <- tl$TL

  plot(net, layout = lay)
}

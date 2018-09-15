################################################
################################################
##### Random Models

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



################################################
################################################
##### Cascade Models


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



################################################
################################################
##### Niche Models


#' Niche Model Food Web
#'
#' @param S Number of species in the community.
#' @param C The connectance, or fraction of realized links in the food web.
#'
#' @return An adjacency matrix for a niche model food web.
#' @export
#'
#' @section Reference:
#' Williams, R. J., and N. D. Martinez. 2000. Simple rules yield complex food webs. Nature 404:180–183.
#'
#' @examples
#' niche(20, .1)
niche <- function(S, C){
  cond <- FALSE
  while(!cond){
    n.i <- sort(runif(S), decreasing = F)
    r.i <- rbeta(S,1,((1/(2*C))-1))*n.i
    c.i <- runif(S, r.i/2, n.i)

    a <- matrix(0, nrow = S, ncol = S)

    for(i in 2:S){
      for(j in 1:S){
        if(n.i[j] > (c.i[i] - (.5 * r.i[i])) & n.i[j] < (c.i[i] + .5 * r.i[i])){
          a[j, i] <- 1
        }
      }
    }

    cond <- igraph::is.connected(igraph::graph.adjacency(a))
  }

  return(a)
}


#' Probabilistic Niche Model Food Web
#'
#' @param S Number of species in the community.
#' @param C The connectance, or fraction of realized links in the food web.
#' @param a The probability that i eats j, when j is exactly on i’s feeding optimum. Default is 0.999, the vaulue used in the original paper.
#'
#' @return An adjacency matrix for a probabilistic niche model food web.
#' @export
#'
#' @section Reference:
#' Williams, R. J., A. Anandanadesan, and D. Purves. 2010. The probabilistic niche model reveals the niche structure and role of body size in a complex food web. PLoS ONE 5.
#'
#' @examples
#' probabilistic_niche(20, .1, .99)
probabilistic_niche <- function(S, C, a = 0.999){
  cond <- FALSE
  while(!cond){
    n.i <- sort(runif(S), decreasing = F)
    r.i <- rbeta(S,1,((1/(2*C))-1))*n.i
    c.i <- runif(S, r.i/2, n.i)

    m <- matrix(0, nrow = S, ncol = S)

    for(i in 2:S){
      for(j in 1:S){
        m[j, i] <- rbinom(1, 1, a*exp(-((n.i[j]-c.i[i])/(r.i[i]/2))^2))
      }
    }

    cond <- igraph::is.connected(igraph::graph.adjacency(m))
  }

  return(m)
}

#' Minimum Potential Niche Model Food Web
#'
#' @param S Number of species in the community.
#' @param C The connectance, or fraction of realized links in the food web.
#'
#' @return An adjacency matrix for a minimum potential niche model food web.
#' @export
#'
#' @section Reference:
#' Allesina, S., D. Alonso, and M. Pascual. 2008. A general model for food web structure. Science (New York, N.Y.) 320:658–61.
#'
#' @examples
min_pot_niche <- function(S, C, f){
  n.i <- sort(runif(S), decreasing = F)
  r.i <- rbeta(S,1,((1/(2*C))-1))*n.i
  c.i <- runif(S, r.i/2, n.i)

  a <- matrix(0, nrow = S, ncol = S)

  for(i in 1:S){
    for(j in 1:S){
      if(n.i[j] > (c.i[i] - (.5 * r.i[i])) & n.i[j] < (c.i[i] + .5 * r.i[i])){
        a[j, i] <- 1
      }
    }
  }
  a <- apply(a, 2, function(x){
    if(sum(x) <= 2){
      return(x)
    }else{
      nlink <- length(x[x == 1])
      x[x == 1][-c(1,nlink)] <- x[x == 1][-c(1,nlink)]*rbinom(nlink-2, 1, prob = (1-f))
      return(x)
    }
  })
  return(a)

}

################################################
################################################
##### Nested Hierarchy Models


#' Nested Hierarchy Food Web Model
#'
#' @param S Number of species in the community.
#' @param C The connectance, or fraction of realized links in the food web.
#'
#' @return An adjacency matrix for a nested hierarchy food web
#' @export
#'
#' @section Reference:
#' Cattin, M., L. Bersier, C. Banašek-Richter, R. Baltensperger, and J.P. Gabriel. 2004. Phylogenetic constraints and adaptation explain food-web structure. Nature 427:835–839.
#' @examples
#' nested_hierarchy(20, .1)
nested_hierarchy <- function(S, C){
  n.i <- (sort(runif(S)))
  l.i <- rbeta(S,1,((1/(2*C))-1))*n.i
  L <- S*(S-1)*C
  l.i <- round((l.i/sum(l.i))*L)
  l.i <- ifelse(l.i > (S-1), S-1, l.i)
  consumers <- which(l.i > 0)
  l.i[1] < -0
  prey <- list()
  for(i in consumers){
    prey[[i]] <- sample(1:(i), 1)
    if(l.i[i] == 1){next}else{
      x <- 2
      while(x <= l.i[i]){
        grp <- sapply(prey[1:(i-1)], function(x) prey[[i]][1] %in% x)
        if(any(grp)){
          p1 <- unlist(prey[grp])
          if(length(p1) >= l.i[i]){prey[[i]] <- c(prey[[i]],sample(p1, l.i[i]-1));x <- l.i[i]+1}else{
            prey[[i]] <- c(prey[[i]],p1)
            x <- length(prey[[i]])+1
          }
          p2 <- (1:(i-1))[-prey[[i]]]
          if(length(p2) > (l.i[i] - length(prey[[i]]))){
            prey[[i]] <- c(prey[[i]], sample(p2, l.i[[i]] - length(prey[[i]])))
            x <- l.i[i] + 1
          }else{
            prey[[i]] <- c(prey[[i]], p2)
            x <- length(prey[[i]])
          }
          if(length(prey[[i]]) < l.i[i]){
            prey[[i]] <- c(prey[[i]], sample(i:length(n.i), (l.i[i]-length(prey[[i]]))))
            x <- length(prey[[i]])+1
          }
        }else{
          p1 <- (1:(i-1))[-prey[[i]]]
          prey[[i]] <- c(prey[[i]], sample(p1), (l.i[i]-length(prey[[i]])))
          x <- length(prey[[i]])+1
        }
      }
    }
  }

  m <- matrix(0, S, S)
  for(i in 1:nrow(m)){
    if(is.null(prey[[i]])){next}
    m[i, prey[[i]]] <- 1
  }

  return(m[rev(order(n.i)),rev(order(n.i))])
}


################################################
################################################
##### Allometric Food Web

#' Allometric Food Web
#'
#' @param S.plant Number of plant species in the community
#' @param S.animal Number of consumer species in the community
#' @param gamma Sets the width of the curve for prey selection. Defaults to 2.
#' @param Ropt Optimal consumer-resource body mass ratio. Defaults to 100.
#' @param thres Minimum feeding efficiency allowed for interactions to occur. Defaults to 0.01.
#'
#' @return
#' @export
#'
#'#' @section Reference:
#' Schneider FD, U Brose, BC Rall, C Guill. 2016. Animal diversity and ecosystem functioning in dynamic food webs. Nature Communications 7:12718.DOI: 10.1038/ncomms12718
#'
#' @examples
allometric <- function(S.plant, S.animal, gamma = 2, Ropt = 100, thres = 0.01){
  mat <- matrix(0, nrow = (S.plant + S.animal), ncol = (S.plant + S.animal))

  log10mass_plant <- 10^runif(S.plant, 0, 6)
  log10mass_animal <- 10^runif(S.animal, 2, 12)

  consumer <- log10mass_animal
  resource <- c(log10mass_plant, log10mass_animal)

  for(i in 1:length(consumer)){
    mat[,(S.plant + i)] <- (consumer[i]/(resource * Ropt) *
                  exp(1 - (consumer[i]/(resource * Ropt))))^2
  }
  mat[mat < thres] <- 0

  return(mat)
}


### Run to create man pages for new/altered functions
###
#
# devtools::document()


cascade_genetic <- function(A){
  S <- nrow(A)
  L <- sum(A)
  v <- lapply(1:10000, function(x) runif(S, 0, 1))
  p.v <- lapply(v, order)

  fits <- sapply(1:length(p.v), function(x){
    m.v <- matrix(0, S, S)
    for(i in 1:S){
      for(j in 1:S){
        if(p.v[[x]][j] == i){m.v[i,j] <- 1}
      }
    }
    A.v <- m.v %*% A %*% t(m.v)
    fitA <- 1/sum(A.v[lower.tri(A.v)])
    return(fitA)
  })
  f2 <- max(fits)
  w <- 1
  while(w < 100){
    nv <- do.call(cbind, v[which(order(fits, decreasing = TRUE) %in% 1:2)])
    nvs <- lapply(1:10000, function(x) apply(nv, 1, function(x) x[sample(c(1,2), 1, replace = T)]))
    p.v <- lapply(nvs, order)

    fits <- sapply(1:length(p.v), function(x){
      m.v <- matrix(0, S, S)
      for(i in 1:S){
        for(j in 1:S){
          if(p.v[[x]][j] == i){m.v[i,j] <- 1}
        }
      }
      A.v <- m.v %*% A %*% t(m.v)
      fitA <- 1/sum(A.v[lower.tri(A.v)])
      return(fitA)
    })
    f2[w] <- max(fits)
    print(f2[w])
    w <- w+1
  }
}

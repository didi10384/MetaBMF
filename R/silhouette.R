ang_dist1 <- function(x){
  distance <- 1-tcrossprod(x)
  diag(distance) <- 0
  as.vector(rowSums(distance))
}

ang_dist2 <- function(x,y){
  distance <- 1-tcrossprod(x,y)
  sum(distance)
}

silhouette <- function(dmat_n,pmat,k)
{
  cl <- cluster(dmat_n,pmat,k)
  clus <- cl$clus
  nbmat <- cl$nbmat
  ts.i <- rep(0,nrow(dmat_n))
  for (j in 1:k){
    Nj <- sum(iC <- clus == j)
    if(Nj != 1){
      dmat_j <- dmat_n[iC,]
      nbmat_j <- nbmat[iC]
      a.i <- ang_dist1(dmat_j)/(Nj - 1)
      b.i <- rep(0,Nj)
      for(i in 1:Nj)
      {
        Ni <- sum(iC1 <- clus == nbmat_j[i])
        dmat_i <- dmat_n[iC1,]
        b.i[i] <- ang_dist2(dmat_j[i,],dmat_i)/Ni
      }
      s.i <- ifelse(a.i != b.i, (b.i - a.i) / pmax(b.i, a.i), 0)
    }
    else{
    s.i <- 0
    }
    ts.i[iC] <- s.i
  }
  mean(ts.i)
}

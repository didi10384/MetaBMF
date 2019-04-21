cluster <- function(dmat_n,pmat,k){
  pmat <- pmat[1:k,]
  bmat <- matrix(0, nrow(dmat_n), nrow(pmat))
  nbmat <- rep(0,nrow(dmat_n))
  ang <- tcrossprod(dmat_n,pmat)
  for(i in 1:nrow(dmat_n)){
    bmat[i,which.max(ang[i,])] <- 1
    nbmat[i] <- which(ang[i,] == sort(ang[i,])[nrow(pmat)-1])
  }
  clus <- which(t(bmat)!=0,arr.ind = T)[,1]
  list(clus=clus,nbmat=nbmat)
}

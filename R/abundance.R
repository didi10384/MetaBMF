abundance <- function(pmat,dmat,clus,clen,nvec){
  norm <- apply(pmat, 1, sum)
  pmat_n <- t(scale(t(pmat), center=FALSE, scale=norm))
  treads <- apply(dmat,1,sum)
  sp_len <- numeric(nrow(pmat))
  sp_tr  <- numeric(nrow(pmat))
  for(i in 1:nrow(pmat)){
    iC <- clus == i
    sp_len[i] <- sum(clen[iC])
    sp_tr[i] <- sum(treads[iC])
  }
  diag(10^6/nvec) %*% t(pmat_n) %*% diag(sp_tr/sp_len * 1000)
}

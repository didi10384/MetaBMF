seeds <- function(dmat,k.max){
  id <- which.max(apply(dmat,1,function(x){sum(x^2)}))
  norm <- apply(dmat, 1, function(x){sqrt(sum(x^2))})
  dmat_n <- t(scale(t(dmat), center=FALSE, scale=norm))
  pmat <- t(dmat_n[id, ])
  pmat1 <- t(dmat[id, ])
  dmat_r <- dmat_n
  dmat_r1 <- dmat
  for(s in 1:(k.max-1)){
    dmat_r <- dmat_r[-id,]
    dmat_r1 <- dmat_r1[-id,]
    angle <- tcrossprod(dmat_r,pmat)
    id <- which.min(apply(angle,1,max))
    pmat <- rbind(pmat,dmat_r[id, ])
    pmat1 <- rbind(pmat1,dmat_r1[id, ])
  }
  list(pmat=pmat,dmat_n=dmat_n,pmat1=pmat1)
}

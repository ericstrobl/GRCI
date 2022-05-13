UpdateCovMat <- function(U,Cov,root){
  
  for (j in setdiff(seq_len(ncol(Cov)),root)){
    Cov[U,j] = (Cov[U,j] - Cov[U,root] * Cov[root,j])/
      (sqrt(1-Cov[U,root]^2)*sqrt(1-Cov[j,root]^2)) 
  }
  
  return(Cov)
  
}
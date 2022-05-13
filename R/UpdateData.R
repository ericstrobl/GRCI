UpdateData <- function(X,U,Cov,root){
  
  for (j in setdiff(U,root)){
    X[,j] = (X[,j] - Cov[j,root] / Cov[root,root] * X[,root,drop=FALSE])/
      sqrt(1-Cov[j,root]^2)
  }
  
  return(X)
  
}
DirectLiNGAM_fast_Y <- function(X,Y,alpha=0.2){
  X = as.matrix(X)
  
  K = c() #1 
  X = normalizeData(X) #3
  Cov = cov(X) #4
  
  U = 1:ncol(X) #2
  
  repeat{ #11
    
    U = U[ttest_fast(X[,U,drop=FALSE],Y,alpha=alpha)] ###
    if (length(U)==0){ ###
      break ###
    } ###
    
    root = FindRoot_fast(X,U,Cov) #6
    K = c(K,root) #7
    U = U[-which(U==root)] #8
    
    if (length(U)==0){ ###
      break ###
    } ###
    X = UpdateData(X,U,Cov,root) #9
    Cov = UpdateCovMat(U,Cov,root) #10
  }
  
  return(list(K=K,X=X[,K,drop=FALSE])) #output
}
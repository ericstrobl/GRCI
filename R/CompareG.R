CompareG <- function(X,i,js,w=NULL){
  
  if (identical(i,js)){  ####
    return(0) ####
  }
  
  rXY = PartialOut(X[,js],X[,i],w[js,,drop=FALSE])
  
  score = c()
  for (k in 1:length(js)){
    score = c(score, one_entropy(X[,js[k]],rXY))
  }
  score = max(score)
  
  return(score) #13  ####
}

one_entropy <- function(X,rXY){
  
  n = length(rXY)
  k = 10
  
  distf = nn2(cbind(X,rXY),k=k+1)$nn.dists[,k+1]
  
  Hx = nn2(as.matrix(X),k=n,searchtype="radius",radius = distf)$nn.idx
  dHx = digamma(rowSums(Hx>0))
  Hx = -mean(dHx)
  
  Hrxy = nn2(as.matrix(rXY),k=n,searchtype="radius",radius = distf)$nn.idx
  dHrxy = digamma(rowSums(Hrxy>0))
  Hrxy = -mean(dHrxy)
  
  return(Hx+Hrxy)
  
}

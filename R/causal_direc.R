causal_direc <- function(X,Y){
  
  rXY = PartialOut_v(X,Y) #9
  rYX = PartialOut_v(Y,X) #10
  
  rXY = rXY / sd(rXY) #11
  rYX = rYX / sd(rYX) #12
  
  X = normalizeData(X)
  Y = normalizeData(Y)
  
  score = diff_entropy(X,Y,rXY,rYX)

  if (score>0){
    return(1)
  } else{
    return(-1)
  }
  
}


diff_entropy <- function(X,Y,rXY,rYX){
  
  n = length(X)
  k = 10
  
  distf = nn2(cbind(X,rXY),k=k+1)$nn.dists[,k+1]
  
  Hx = nn2(as.matrix(X),k=n,searchtype="radius",radius = distf)$nn.idx
  dHx = digamma(rowSums(Hx>0))
  
  Hrxy = nn2(as.matrix(rXY),k=n,searchtype="radius",radius = distf)$nn.idx
  dHrxy = digamma(rowSums(Hrxy>0))
  
  score1 = -mean(dHx+dHrxy)
  # /sd(dHx+dHrxy)
  
  distf = nn2(cbind(Y,rYX),k=k+1)$nn.dists[,k+1]
  
  Hy = nn2(as.matrix(Y),k=n,searchtype="radius",radius = distf)$nn.idx
  dHy = digamma(rowSums(Hy>0))
  
  Hryx = nn2(as.matrix(rYX),k=n,searchtype="radius",radius = distf)$nn.idx
  dHryx = digamma(rowSums(Hryx>0))
  
  score2 = -mean(dHy+dHryx)

  return(-score1+score2)

}

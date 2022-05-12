LiNGAM_direc <- function(X,Y){
  
  rXY = lm.fit(cbind(1,X),Y)$residuals #9
  rYX = lm.fit(cbind(1,Y),X)$residuals #10
  
  rXY = rXY / sd(rXY) #11
  rYX = rYX / sd(rYX) #12
  
  X = normalizeData(X)
  Y = normalizeData(Y)
  
  score = LRT2(X,Y,rXY,rYX)
  if (score>0){
    return(1)
  } else{
    return(-1)
  }
  
}


LRT2 <- function(X,Y,rXY,rYX){
  
  return(-Hu(X)-Hu(rXY)+Hu(Y)+Hu(rYX))
}

Hu <- function(u){
  
  k1 = 79.047
  k2 = 7.4129
  beta = 0.37457
  
  H = 0.5*(1+log(2*pi))-k1*(mean(log(cosh(u)))-beta)^2 - k2*mean(u*exp(-(u^2)/2))^2
  
  return(H)
}

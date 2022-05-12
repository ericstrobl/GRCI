Compare2 <- function(X,i,j,Cov){
  
  if (i==j){  ####
    return(0) ####
  }
  
  rij = X[,i] - Cov[i,j]/Cov[j,j] * X[,j] #9
  rji = X[,j] - Cov[j,i]/Cov[i,i] * X[,i] #10
  
  rsd = sqrt(1-Cov[i,j]^2)
  rij = rij / rsd #11
  rji = rji / rsd #12
  
  score = LRT(X[,i],X[,j],rij,rji)
  
  return(score) #13  ####
}

LRT <- function(xi,xj,rij,rji){
  
  return(Hu(xj)+Hu(rij)-Hu(xi)-Hu(rji))
}

Hu <- function(u){
  
  k1 = 79.047
  k2 = 7.4129
  beta = 0.37457
  
  H = 0.5*(1+log(2*pi))-k1*(mean(log(cosh(u)))-beta)^2 - k2*mean(u*exp(-(u^2)/2))^2
  
  return(H)
}
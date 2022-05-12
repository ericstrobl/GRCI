FindRoot_fast <- function(X,U,Cov){
  
  r = length(U) #4
  
  if (r==1){ #1
    return(U) #2
  }
  
  O = rep(1,r)
  S = rep(0,r)
  M = matrix(TRUE,r,r)
  
  repeat{
    I = (S == min(S))
    
    if ( sum(O[I] > r) > 0){
      break
    }
    
    for (i in which(I)){
      if (M[i,O[i]]){
        score = Compare2(X,U[i],U[O[i]],Cov)
        S[i] = S[i] + min(0,score)^2
        S[O[i]] = S[O[i]] + min(0,-score)^2
        
        M[i,O[i]]=FALSE
        M[O[i],i]=FALSE
      }
      
      O[i] = O[i]+1
    }
  }
  
  root = U[S==min(S)][1]
  
  return(root) #output
  
}
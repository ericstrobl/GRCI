sample_DAG_HNM <- function(nsamps, DAG){
  
  G = DAG$graph
  r = nrow(G)
  
  Y = DAG$Y
  
  err=matrix(0,nsamps,r)
  for (i in 1:r){
    if (DAG$errors[i] == 1){
      err[,i]=matrix(rt(nsamps,df=5),nsamps)
    } else if (DAG$errors[i]==2){
      err[,i]=matrix(runif(nsamps,-1,1),nsamps)
    } else if (DAG$errors[i]==3){
      err[,i]=matrix(rchisq(nsamps,df=3)-3,nsamps)
    }
  }
  
  
  err[,Y]=0 #error for diagnosis is zero
  data=normalizeData(err)
  
  done=which(colSums(G)==0) # variables without parents
  stop=0;
  while (stop==0){
    for (s in done){
      ch=which(G[s,]==1) 
      for (c in ch){
        if (c %in% done){
          next
        }
        pa=which(G[,c]==1) 
        
        h=intersect(pa,done)
        
        if (setequal(h,pa)){ # if all parents already done
          
          M = (data[,h,drop=FALSE]%*%DAG$weightsM[h,c,drop=FALSE])
          M = (M - mean(M))/sd(M)
          if (DAG$functionsM[c]==2){ # multiplicative effect
            errc = err[,c]*(sqrt(M^2+1)) # shift to allow to linear correlation
          } else if(DAG$functionsM[c]==3){
            errc = err[,c]*(M*pnorm(M) + 1)
          } else if(DAG$functionsM[c]==4){
            errc = err[,c]*(1/(1+exp(-M)) + 1)
          }
          
          # errc = err[,c]
  
          A = (data[,h,drop=FALSE]%*%DAG$weightsA[h,c,drop=FALSE])
          A = (A - mean(A))/sd(A)
          if (DAG$functionsA[c]==2){  # additive effect
            A = sqrt(A^2+1) - 1
            # A = log(1+exp(A))
            # A = sinh(A)
          } else if(DAG$functionsA[c]==3){
            A = A*pnorm(A)
          } else if(DAG$functionsA[c]==4){
            A = 1/(1+exp(-A))
          }
          data[,c]=A+errc
   
          done=unique(c(done, c))
        }
      }
    }
    
    if (length(done) == r){
      stop=1;
    }
  }
  
  Y0 = (data[,Y]-mean(data[,Y]))/sd(data[,Y])
  pY = logistic(Y0)
  for (i in 1:nsamps){
    data[i,Y] = rbinom(n=1,size=1,prob=pY[i]) 
  }
  
  return(list(data=data,E=err,Y0=Y0))
}


logistic <- function(X){
  
  return(1/(1+exp(-X)))
}


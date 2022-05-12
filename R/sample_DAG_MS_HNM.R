sample_DAG_MS_HNM <- function(X,Y,G,E,w,outm=NULL){
  
  # new X to be filled in
  Xn = matrix(0,nrow(X),ncol(X)+1)
  Xn[,1:ncol(X)] = E
  
  if (is.null(outm)){
    outm = vector(mode = "list", length = (ncol(X)+1))
    mstart = TRUE
  } else{
    mstart = FALSE
  }
  
  X = as.matrix(X)
  n = nrow(X)
  Xa = normalize01(rbind(X,Xn[,1:ncol(X)])) ##
  X = Xa[1:n,,drop=FALSE] ##
  Xn = Xa[(n+1):nrow(Xa),,drop=FALSE] ##
  Xn = cbind(Xn,0)

  means = Xn
  done=which(colSums(G)==0)
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
          if (length(h)>0){
            if (c != (ncol(X)+1)){
              if (mstart){
                outm[[c]] = PredictOut(X[,h,drop=FALSE],X[,c,drop=FALSE],Xn[,h,drop=FALSE],w[h,,drop=FALSE])
                SD = outm[[c]]$SD
                Mean = outm[[c]]$Mean
              } else{
                SD = make_features(Xn[,h,drop=FALSE],outm[[c]]$k2,w[h,,drop=FALSE]) %*% outm[[c]]$modSD
                Mean = make_features(Xn[,h,drop=FALSE],outm[[c]]$k1,w[h,,drop=FALSE]) %*% outm[[c]]$modMean
              }
              Xn[,c] = SD*E[,c,drop=FALSE] + Mean
            } else{
              if (mstart){
                outm[[c]] = PredictOut(X[,h,drop=FALSE],as.matrix(Y),Xn[,h,drop=FALSE],w[h,,drop=FALSE],binary_tar = TRUE)
                Mean = outm[[c]]$Mean
              } else{
                Mean = make_features(Xn[,h,drop=FALSE],outm[[c]]$k,w[h,,drop=FALSE]) %*% outm[[c]]$modMean
              }
              Xn[,c] = Mean
            }
          } else{
            if (c != (ncol(X)+1)){
              Xn[,c] = E[,c]
            } else{
              Xn[,c] = Y
            }
          }
          
          done=unique(c(done, c))
        }
      }
    }
    
    if (length(done) == (ncol(X)+1)){
      stop=1;
    }
  }
  
  
  return(list(X = Xn,E=E,outm=outm))
}

PredictOut <- function(X,Y,X_test,w=NULL,binary_tar = FALSE){
  
  Y = as.matrix(Y)
  X = as.matrix(X)
  
  if (binary_tar){
    
    Res = spline_regressionR_target_test(X,Y,X_test,w=w)
    return(list(SD=0, Mean = Res$Yhat, modMean = Res$mod, k = Res$k))
    
  } else{
    ## MEAN
    Res = spline_regressionR_test(X,Y,X_test,w=w)
    
    ## VARIANCE
    SD = spline_regressionR_test(X,abs(Y-Res$Yhat_cv),X_test,w = w, pos=TRUE)
    
    return(list(SD = SD$Yhat, Mean = Res$Yhat, modMean = Res$mod, modSD = SD$mod,
                k1 = Res$k, k2 = SD$k))
    
  }
}

make_features <- function(X,k,w){
  n = nrow(as.matrix(X))
  k_num = unique(round(seq(2,min(200,sqrt(n/10)),length.out=10)))
  if (k==0){
    Xd = matrix(1,n,1) # constant
  } else if (k==1){ 
    Xd = cbind(X,1) # linear
  } else{
    kX = seq(0,1,length.out=k_num[k])
    kX = t(replicate(n, kX))
    Xd = as.matrix(pmax(X %*% w[,1:k_num[k]]-kX,0)) # spline
    Xd = Xd[,-c(1,k_num[k])] # remove first and last kX
    Xd = cbind(Xd,1)  # add in constant term for last kX
    Xd = cbind(X,Xd) # add in linear terms for first kX
  }
  
  return(Xd)
}


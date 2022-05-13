spline_regressionR_binary <- function(X,Y,w=NULL,penalty = 1){
  
  #higher penalty = more multivariate-ness
  
  X = normalize01(as.matrix(X))
  Y = as.matrix(Y)
  Y = cbind(Y, (Y-1)*-1)
  
  n = nrow(X)
  
  k_num = unique(round(seq(2,min(200,sqrt(n/10)),length.out=10)))
  
  if (is.null(w)){
    w = rdirichlet(200,rep(1,ncol(X))*penalty) # random projections
    # w = normalize_rep(w) # weight each feature equally
    w = t(w)
  } else{
    w = w/colSums(w)
  }
  
  err=Inf
  for(k in 0:length(k_num)){
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
    
    H = Xd %*% spdinv(t(Xd) %*% Xd + diag(ncol(Xd))*1E-10) %*% t(Xd)
    Yhat_t = H %*% Y
    
    Pre_t = Y-(Y-Yhat_t)/(1-diag(H))
    
    Yhat_t = pmax(Yhat_t,0) + 1E-10###
    Yhat_t = Yhat_t/rowSums(Yhat_t)
    
    Pre_t = pmax(Pre_t,0) + 1E-10###
    Pre_t = Pre_t/rowSums(Pre_t)
    
    err_t = mean( (Y-Pre_t)^2 )
    
    if (err_t<err){
      err = err_t
      Yhat = Yhat_t
      Yhat_cv = Pre_t
    }
    
  }
  
  return(list(Yhat = Yhat[,1,drop=FALSE], Yhat_cv = Yhat_cv[,1,drop=FALSE], w = w))
  
  
}
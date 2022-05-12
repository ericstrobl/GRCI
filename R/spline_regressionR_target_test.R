spline_regressionR_target_test <- function(X,Y,X_test,w=NULL,penalty = 1){
  # spline regression, assumes target is non-negative
  
  X = as.matrix(X)
  n = nrow(X)
  
  X_test = as.matrix(X_test)
  
  Y = as.matrix(Y)
  sY = sd(Y)
  Y = t(t(Y)/sY)
  
  k_num = unique(round(seq(2,min(200,sqrt(n/10)),length.out=10)))
  
  if (is.null(w)){
    require(DirichletReg)
    w = rdirichlet(200,rep(1,ncol(X))*penalty) # random projections
    # w = normalize_rep(w) # weight each feature equally
    w = t(w)
  } else{
    w = w/colSums(w)
  }
  
  nt = nrow(X_test)
  
  err=Inf
  for(k in 0:length(k_num)){
    if (k==0){
      Xd = matrix(1,n,1) # constant
      
      Xe = matrix(1,nt,1)
    } else if (k==1){ 
      Xd = cbind(X,1) # linear
      
      Xe = cbind(X_test,1)
    } else{
      kX = seq(0,1,length.out=k_num[k])
      kX = t(replicate(n, kX))
      Xd = as.matrix(pmax(X %*% w[,1:k_num[k]]-kX,0)) # spline
      Xd = Xd[,-c(1,k_num[k])] # remove first and last kX
      Xd = cbind(Xd,1)  # add in constant term for last kX
      Xd = cbind(X,Xd) # add in linear terms for first kX
      
      kX = seq(0,1,length.out=k_num[k])
      kX = t(replicate(nt, kX))
      Xe = as.matrix(pmax(X_test %*% w[,1:k_num[k]]-kX,0)) # spline
      Xe = Xe[,-c(1,k_num[k])] # remove first and last kX
      Xe = cbind(Xe,1)  # add in constant term for last kX
      Xe = cbind(X_test,Xe) # add in linear terms for first kX
    }
    
    Hi = spdinv(t(Xd) %*% Xd + diag(ncol(Xd))*1E-10) %*% t(Xd)
    Htr = Xd %*% Hi ##
    Hte = Xe %*% Hi ##
    
    Yhat_t = Htr %*% Y ##
    Yhat_t = pmax(Yhat_t,0) ## non-negative
    
    Pre_t = Y-(Y-Yhat_t)/(1-diag(Htr)) ##
    Pre_t = pmax(Pre_t,0) ## non-negative
    
    err_t = colMeans( (Y-Pre_t)^2 ) ##
    
    # print(err_t)
    
    if (err_t<err){
      err = err_t
      kf = k
      mod = Hi %*% Y
      Yhat = Hte %*% Y ##
      Yhat_cv = Pre_t ##
    }
    
  }
  
  return(list(Yhat = Yhat*sY, Yhat_cv = Yhat_cv*sY, w = w, mod = mod, k = kf))
  
  
  
}
spline_regressionR_test <- function(X,Y,X_test,w=NULL,pos=FALSE, penalty = 1){
  
  X = as.matrix(X)
  n = nrow(X)
  
  X_test = as.matrix(X_test) ##
  
  Y = as.matrix(Y)
  c = ncol(Y)
  
  sY = sd(Y)
  Y = t(t(Y)/sY)
  
  n = nrow(X)
  
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
    
    Yhat_t = Htr %*% Y
    if (pos){
      Yhat_t = pmax(Yhat_t,5E-2) # SD must be >0
    }
    Res_t = (Y-Yhat_t)/(1-diag(Htr))
    err_t = colMeans( Res_t^2 )
    
    # print(err_t)
    
    if (err_t<err){
      err = err_t
      Yhat = Hte %*% Y
      mod = Hi %*% Y
      kf = k
      if (pos){
        Yhat_cv = pmax(Y-Res_t,5E-2)
      } else{
        Yhat_cv = Y-Res_t
      }
    }
    
  }
  
  # print(ncol(Hte))
  
  return(list(Yhat = Yhat*sY, Yhat_cv = Yhat_cv*sY, w = w, mod = mod, k = kf))
  
  
}
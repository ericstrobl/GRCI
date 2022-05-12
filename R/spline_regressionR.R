spline_regressionR <- function(X,Y,w=NULL,pos=FALSE, penalty = 1){
  
  #higher penalty = more multivariate-ness
  
  X = normalize01(as.matrix(X))
  Y = as.matrix(Y)
  c = ncol(Y)
  
  sY = sd(Y)
  Y = t(t(Y)/sY)
  
  n = nrow(X)
  
  k_num = unique(round(seq(2,min(200,sqrt(n/10)),length.out=10))) ###
  
  if (is.null(w)){
    w = rdirichlet(200,rep(1,ncol(X))*penalty) # random projections
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
    
    H = Xd %*% spdinv(t(Xd) %*% Xd + diag(ncol(Xd))*1E-8) %*% t(Xd) #### was 1E-10
    Yhat_t = H %*% Y
    if (pos){
      Yhat_t = pmax(Yhat_t,5E-2) # SD must be >0
    }
    Res_t = (Y-Yhat_t)/(1-diag(H))
    err_t = colMeans( Res_t^2 )
    
    # print(err_t)
    
    if (err_t<err){
      kf = k
      err = err_t
      Yhat = Yhat_t
      if (pos){
        Yhat_cv = pmax(Y-Res_t,5E-2)
      } else{
        Yhat_cv = Y-Res_t
      }
    }
    
  }
  
  # print(kf)
  
  return(list(Yhat = Yhat*sY, Yhat_cv = Yhat_cv*sY, w = w))
  
  
}

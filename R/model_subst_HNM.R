model_subst_HNM <- function(X,Y,G,E){
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  # plot(as(G,"graphNEL"))
  
  i0 = which(Y==0)
  i1 = which(Y==1)
  n0 = length(i0)
  E0 = E[i0,]
  
  require(DirichletReg)
  
  w = rdirichlet(200,rep(1,ncol(X))) # random projections
  w = t(w)
  
  c = ncol(X)
  pa_Y = which(G[,c+1]>0)
  if (length(pa_Y)>0){
    Yp = spline_regressionR_binary(X[,pa_Y,drop=FALSE],Y,w=w[pa_Y,,drop=FALSE])$Yhat_cv
  } else{
    return(list(order=NULL,scores=0))
  }
  
  Anc = c()
  score = c()
  outS = sample_DAG_MS_HNM(X[i0,],Yp[i0],G,E0,w)
  L0 = mean(outS$X[,c+1]) # score for healthy
  for (i in 1:c){ # compute marginal contributions
    if (isAnc(G,i,c+1)){
      Anc = c(Anc, i)
      
      En = E0
      En[,i]= sample(E[i1,i],n0,replace=TRUE) # bootstrap substitute the errors for diseased
      
      Ln = sample_DAG_MS_HNM(X[i0,],Yp[i0],G,En,w,outS$outm)$X # sample new DAG
      
      score = c(score,mean(Ln[,c+1]) - L0)
    }
  }
  
  return(list(order=Anc,scores=score)) # scores aligned with Anc
  
  
}
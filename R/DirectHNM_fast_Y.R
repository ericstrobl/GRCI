DirectHNM_fast_Y <- function(X,Y,G,alpha=0.2){
  X = as.matrix(X)
  # E = X
  
  K = c()
  S = rep(-Inf,ncol(X))
  U = 1:ncol(X)
  update = U
  
  penalty = 1
  w = rdirichlet(200,rep(1,ncol(X))*penalty) # random projections
  # w = normalize_rep(w)
  w = t(w)
  
  repeat{ 
    s_out = FindSink(X,U,S,G,update,w)
    sink = s_out$sink
    S = s_out$S
    
    K = c(K,sink)
    U = U[-which(U==sink)]
    
    if (length(U)==0){ ###
      break ###
    } ###
    
    if (sum(G[sink,])){
      X[,sink] = PartialOut(X[,G[sink,]],X[,sink],w[G[sink,],,drop=FALSE]) # partial out neighbors from sink
    }
    update = intersect(U,which(G[sink,]))
    G[sink,] = FALSE; G[,sink]=FALSE # remove node
  
  }
  
  return(list(K=K,X=X[,K,drop=FALSE])) #output
}
FindSink <- function(X,U,S,G,update,w=NULL){
  r = length(update)
  
  ## BASELINE SCORES
  for (i in seq_len(r)){
    V = which(G[update[i],])
    if (length(V)==0){ # if no neighbors, then break because score is -Inf
      next
    }

    S[update[i]] = CompareG(X,update[i],V,w) # baseline

  }
  
  sink = U[S[U]==min(S[U])][1] 
  
  return(list(sink=sink,S=S)) #output
  
}

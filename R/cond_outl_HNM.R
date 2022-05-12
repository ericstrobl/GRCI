cond_outl_HNM <- function(X,Y,G){
  
  Anc = c()
  c = ncol(X)
  CO = matrix(0,nrow(X),c)
  for (i in 1:c){
    if (isAnc(G,i,c+1)){
      Anc = c(Anc, i)
      
      pa = which(G[,i]>0)
      if (length(pa)>0){
        CO[,i] = abs(PartialOut(X[,pa,drop=FALSE],X[,i,drop=FALSE]))
      } else{
        CO[,i] = abs(X[,i]-mean(X[,i]))/sd(X[,i])
      }
    }
  }
  
  return(list(order=Anc,scores=CO[,Anc,drop=FALSE]))
  
}

isAnc <- function(graph,a,b,visited=rep(FALSE,nrow(graph)))
{
  
  if (a %in% b){
    return(TRUE)
  }
  
  visited[a] = TRUE;
  
  adj = which(graph[a,] & !visited);
  
  out=FALSE;
  for (j in adj){
    out=isAnc(graph, j, b, visited);
    if(out==TRUE){
      break;
    }
  }
  
  return(out)
  
}



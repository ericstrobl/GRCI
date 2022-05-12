FindAnc <- function(graph, tar){
  
  ind = 1:ncol(graph)
  ind = ind[-tar]
  
  Anc = c()
  for (i in ind){
     if (isAnc_fast_LE(graph,i,tar)){
       Anc = c(Anc, i)
     }
  }
  
  return(Anc)
}



isAnc_fast_LE <- function(graph,a,b,visited=rep(FALSE,nrow(graph)))
{
  
  if (a %in% b){
    return(TRUE)
  }
  
  visited[a] = TRUE;
  
  adj = which(graph[a,] & !visited);
  
  out=FALSE;
  for (j in adj){
    out=isAnc_fast_LE(graph, j, b, visited);
    if(out==TRUE){
      break;
    }
  }
  
  return(out)
  
}

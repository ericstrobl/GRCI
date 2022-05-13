generate_DAG_HNM <- function(p,en){
  
  N = p*p - p;
  samplesB = rbinom(N/2,1, en/(p-1) );
  
  DAG = list()
  graph = matrix(0,p,p)
  graph[upper.tri(graph, diag=FALSE)] <- samplesB;
  
  DAG = list()
  DAG$graph = graph

  weightsA = matrix((1*runif(p^2)+0.5)*sample(c(-1,1),p^2,replace=TRUE),p,p)
  # weightsA = matrix(1,p,p)
  DAG$weightsA = weightsA*DAG$graph
  
  weightsM = matrix((1*runif(p^2)+0.5)*sample(c(-1,1),p^2,replace=TRUE),p,p)
  # weightsM = matrix(1,p,p)
  DAG$weightsM = weightsM*DAG$graph
  
  ord = sample(1:p,p,replace=FALSE) # permute order
  DAG$graph = DAG$graph[ord,ord]
  DAG$weights = DAG$weights[ord,ord]
  DAG$weightsA = DAG$weightsA[ord,ord]
  DAG$weightsM = DAG$weightsM[ord,ord]
  
  DAG$functionsA = sample(2:4,p,replace=TRUE)
  DAG$functionsM = sample(2:4,p,replace=TRUE)
  
  DAG$errors = sample(1:3,p,replace=TRUE)
  
  pt = which(rowSums(DAG$graph)==0 & colSums(DAG$graph)>0)
  if (length(pt)>1){
    DAG$Y = sample(which(rowSums(DAG$graph)==0 & colSums(DAG$graph)>0),1) #terminal vertex with at least one parent
  } else{
    DAG$Y = pt
  }
  
  return(DAG)
}
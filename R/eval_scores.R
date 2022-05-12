eval_scores <- function(true_scores,output){
  
  if (length(output$order)==0){
    return(0)
  }
  est_scores = matrix(-Inf,nrow(true_scores),ncol(true_scores))
  est_scores[,output$order]= output$scores
  
  n = nrow(true_scores)
  vars = 1:ncol(true_scores)
  SIMs = c()
  for(i in seq_len(n)){
    ix = which(true_scores[i,]>0)
    oix = order(true_scores[i,],decreasing=TRUE)[1:length(ix)]
    
    true_sort = vars[oix]
    est_sort = vars[order(est_scores[i,],decreasing=TRUE)]
    SIMs = c(SIMs,SIM(true_sort,est_sort,true_scores[i,oix]))
    # SIMs = c(SIMs,SIM(true_sort,est_sort,rep(1,length(true_scores[i,oix]))))
  }
  return(mean(SIMs))
}

SIM <- function(true_sort,est_sort,weights){
  if (sum(weights)==0){
    return(1)
  } else{
    weights = weights/sum(weights)
  }
  d = length(true_sort)
  ni = 0
  for (i in seq_len(d)){
    A = length(intersect(true_sort[1:i],est_sort[1:i]))/i
    ni = ni + weights[i]*A
  }
  return(ni)
}

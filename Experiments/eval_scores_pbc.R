eval_scores_pbc <- function(output,Xo){
  
  est_scores = matrix(-Inf,nrow(Xo),ncol(Xo))
  est_scores[,output$order]= output$scores
  n = nrow(Xo)
  vars = 1:ncol(Xo)
  SIMs = c()
  
  for(i in seq_len(n)){
    if (Xo[i,2]<2){
      oix = c(1,2)
    } else{
      oix = c(2,1)
    }
    
    true_sort = vars[oix]
    est_sort = vars[order(est_scores[i,],decreasing=TRUE)]
    SIMs = c(SIMs,SIM(true_sort,est_sort,c(1,1)))
  }
  return(mean(SIMs))
}

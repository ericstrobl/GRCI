eval_scores_diabetes <- function(output,Xo,id){
  
  est_scores = matrix(-Inf,nrow(Xo),ncol(Xo))
  est_scores[,id[output$order]]= output$scores
  n = nrow(Xo)
  vars = 1:ncol(Xo)
  SIMs = c()
  
  for(i in seq_len(n)){
    oix = order(c(0,Xo[i,2],rep(0,5),Xo[i,8]),rnorm(8),decreasing=TRUE)
    if (Xo[i,2]>=140){
      oix = order(c(0,1,rep(0,6)),rnorm(8),decreasing=TRUE)
    } else{
      oix = order(c(0,Xo[i,2],rep(0,5),Xo[i,8]),rnorm(8),decreasing=TRUE)
    }
    
    true_sort = vars[oix]
    est_sort = vars[order(est_scores[i,],decreasing=TRUE)]
    SIMs = c(SIMs,SIM(true_sort,est_sort,rep(1,length(oix))))
  }
  return(mean(SIMs))
}

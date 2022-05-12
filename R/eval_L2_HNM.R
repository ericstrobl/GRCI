eval_L2_HNM <- function(true_scores,output_E,output_order,true_anc){
  
  est_scores = matrix(0,nrow(true_scores),ncol(true_scores))
  if (!is.null(output_E)){
    output_E = t(output_E) - colMeans(output_E)
    output_E  = t(output_E / apply(output_E,1,function(x) mean(abs(x))))
    est_scores[,output_order]= output_E
  }
  
  # print(colMeans(true_scores[,true_anc,drop=FALSE]))
  # print(est_scores[1:5,true_anc])
  
  return(mean( (true_scores[,true_anc] - est_scores[,true_anc])^2 ))
}
normalizeData_HNM <- function(X){
  
  X = as.matrix(X)
  for (i in seq_len(ncol(X))){
    # if (sd(X[,i]) == 0){
    diff = X[,i]-mean(X[,i])
    mad = mean(abs(diff))
    if (mad==0){
      X[,i] = diff
    } else{
      # X[,i] = (X[,i] - mean(X[,i]))/sd(X[,i])
      X[,i] = diff/mad
    }
  }
  X = as.matrix(X)
  
  return(X)
}
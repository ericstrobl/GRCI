PartialOut <- function(X,Y,w=NULL){
  
  Y = as.matrix(Y)
  X = as.matrix(X)

  n = length(X)

  ## CONDITIONAL MEAN
  Res = spline_regressionR(X,Y,w=w)

  ### CONDITIONAL MEAN ABSOLUTE DEVIATION OF THE MEAN
  MAD = spline_regressionR(X,abs(Y-Res$Yhat_cv),w = w, pos=TRUE)

  Res = (Y - Res$Yhat)/MAD$Yhat

  return(Res)
}

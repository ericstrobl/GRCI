ttest_fast <- function(X,Y,alpha=0.2){
  
  X0 = X[Y==0,,drop=FALSE]
  X1 = X[Y==1,,drop=FALSE]
  
  ms = colMeans(X0) - colMeans(X1)
  
  n0 = nrow(X0)
  n1 = nrow(X1)
  
  var0 = apply(X0,2,var)
  var1 = apply(X1,2,var)
  
  t = ms/ sqrt(var0/n0 + var1/n1)
  
  df = ((var0 + var1)^2) / ((var0^2)/(n0-1) + (var1^2)/(n1-1))
  
  ps = (1-pt(abs(t),df=df))*2
  
  return( which(ps<alpha) )
  
  
  
}
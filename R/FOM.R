FOM <- function(X,Y){
  
  X = normalize01(X)
  Y = normalize01(Y)
  
  # X_p = as.matrix(seq(0,1,length.out=1000))
  nvar = 1
  m1 = mleHetGP(X = X, Z = Y, lower = rep(0.1, nvar), upper = rep(50, nvar),
                         covtype = "Matern5_2")
  m1 = predict(x=X, object = m1)
  for1 = mean(3*(m1$sd2 + m1$nugs)^2)
  
  m2 = mleHetGP(X = Y, Z = X, lower = rep(0.1, nvar), upper = rep(50, nvar),
                covtype = "Matern5_2")
  m2 = predict(x=Y, object = m2)
  back1 = mean(3*(m2$sd2 + m2$nugs)^2)
  
  # print(for1)
  # print(back1)
  
  if (for1<back1){
    return(1)
  } else{
    return(-1)
  }
}
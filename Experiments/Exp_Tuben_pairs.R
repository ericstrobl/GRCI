
library(MASS)
library(RANN)
library(Rfast)
library(hetGP)
library(gptk)

results = list()
skip = c(47, 52, 53, 54, 55, 70, 71, 105, 107)
load("pairs.RData")
for (i in 1:108){
  print(i)
  if (is.element(i,skip)){
    next
  }
  
  X = pairs[[i]]$X[,1,drop=FALSE]
  Y = pairs[[i]]$Y[,1,drop=FALSE]
  results$total = c(results$total, pairs[[i]]$weight)

  #GRCI
  ptm <- proc.time()
  out = causal_direc(X,Y)
  results$GRCI$time = c(results$GRCI$time, (proc.time() - ptm)[3])
  if (out >= 0){
    results$GRCI$scores = c(results$GRCI$scores, pairs[[i]]$weight)
  } else{
    results$GRCI$scores = c(results$GRCI$scores, 0)
  }

  #HEC
  ptm <- proc.time()
  out = HEC(X,Y)
  results$HEC$time = c(results$HEC$time, (proc.time() - ptm)[3])
  if (out >= 0){
    results$HEC$scores = c(results$HEC$scores, pairs[[i]]$weight)
  } else{
    results$HEC$scores = c(results$HEC$scores, 0)
  }

  #FOM
  n = min(nrow(X),1000)
  ptm <- proc.time()
  out = FOM(X[1:n,1],Y[1:n,1])
  results$FOM$time = c(results$FOM$time, (proc.time() - ptm)[3])
  if (out >= 0){
    results$FOM$scores = c(results$FOM$scores, pairs[[i]]$weight)
  } else{
    results$FOM$scores = c(results$FOM$scores, 0)
  }

  n = min(nrow(X),3000)
  ptm <- proc.time()
  out <- ICML(cbind(normalize01(X[1:n,1]),normalize01(Y[1:n,1])), model = train_gp, indtest = indtestHsic, output = FALSE)
  results$ANM$time = c(results$ANM$time, (proc.time() - ptm)[3])
  if (identical(out$Cd, "->")){
    results$ANM$scores = c(results$ANM$scores, pairs[[i]]$weight)
  } else{
    results$ANM$scores = c(results$ANM$scores, 0)
  }

  #Direct LiNGAM
  ptm <- proc.time()
  out = LiNGAM_direc(X,Y)
  results$LiNGAM$time = c(results$LiNGAM$time, (proc.time() - ptm)[3])
  if (out >= 0){
    results$LiNGAM$scores = c(results$LiNGAM$scores, pairs[[i]]$weight)
  } else{
    results$LiNGAM$scores = c(results$LiNGAM$scores, 0)
  }

  save(file="Res_Tuben_pairs2.RData",results)
  
}

# accuracy
plot(1:length(results$GRCI$scores),cumsum(results$GRCI$scores)/cumsum(results$total),type="l",ylim=c(0,1))

# timing
times = matrix(0,reps,7)
for (i in 1:reps){
  times[i,1] = results$GRCI$time
  times[i,4] = results$HEC$time
  times[i,5] = results$FOM$time
  times[i,2] = results$ICML$time
  times[i,3] = results$LiNGAM$time
}

print(colMeans(times))
print(apply(times,2,sd)/sqrt(reps)*1.96)


library(MASS)
library(RANN)
library(Rfast)
library(hetGP)
library(gptk)

nsamps = 1000
results = list()
results$GRCI$time = matrix(0,5,500)
results$GRCIA$time = matrix(0,5,500)
results$HEC$time = matrix(0,5,500)
results$FOM$time = matrix(0,5,500)
results$ANM$time = matrix(0,5,500)
results$LiNGAM$time = matrix(0,5,500)

results$GRCI$scores = matrix(0,5,500)
results$GRCIA$scores = matrix(0,5,500)
results$HEC$scores = matrix(0,5,500)
results$FOM$scores = matrix(0,5,500)
results$ANM$scores = matrix(0,5,500)
results$LiNGAM$scores = matrix(0,5,500)

for (j in 1:5){
  for (i in 1:200){
    print(i)
    
    if (j %in% 1:3){ # non-Gaussian error for X
      xt = sample(1:3,1) 
      if (xt == 1){
        X = runif(nsamps)*2-1
      } else if (xt == 2){
        X = rt(nsamps,df=5)
      } else if (xt == 3){
        X = rchisq(nsamps,df=3)
      }
    } else{
      X = rnorm(nsamps) # Gaussian error for X
    }
    
    if (j %in% 1:3){ # non-Gaussian error for Y
      yt = sample(1:3,1)
      if (yt == 1){
        Y = runif(nsamps)*2-1
      } else if (yt == 2){
        Y = rt(nsamps,df=5)
      } else if (yt == 3){
        Y = rexp(nsamps)
      }
    } else{
      Y = (1/3)*rnorm(nsamps)  # Gaussian error for Y. 0.5 is approximately mean of 1/5 and sqrt(2/5) from page 7 of Tagasovska et al. (2020)
    }
    
    if (j == 1){ # linear case
      Y = Y + X*(0.75*runif(1)+0.25)*sample(c(-1,1),1) #coefficient in U([-1,-0.25] \cup [0.25, 1])
    } else if ((j==3) | (j==5)){ # heteroscedastic cases
      yf1 = sample(1:3,1)
      if (yf1 == 1){
        g = sqrt(X^2+1)
        Y = Y * g
      } else if (yf1 == 2){
        g = X*pnorm(X) + 1
        Y = Y * g
      } else if (yf1 == 3){
        g = 1/(1+exp(-X)) + 1
        Y = Y * g
      }
    }
    
    if (j!=1){ # non-linear cases
      yf = sample(1:3,1)
      if (yf == 1){
        Y = Y + sqrt(X^2+1)-1
      } else if (yf == 2){
        Y = Y + X*pnorm(X)
      } else if (yf == 3){
        Y = Y + 1/(1+exp(-X))
      }
    }
    
    X = normalizeData(X) # no gaming of variances
    Y = normalizeData(Y)
    
    # GRCI
    ptm <- proc.time()
    out = causal_direc(X,Y)
    results$GRCI$time[j,i] = (proc.time() - ptm)[3]
    if (out >= 0){
      results$GRCI$scores[j,i] = 1
    } else{
      results$GRCI$scores[j,i] = 0
    }
    
    # #GRCI, additive noise model
    # ptm <- proc.time()
    # out = causal_direc_ANM(X,Y)
    # results$GRCIA$time[j,i] = (proc.time() - ptm)[3]
    # if (out >= 0){
    #   results$GRCIA$scores[j,i] = 1
    # } else{
    #   results$GRCIA$scores[j,i] = 0
    # }
    
    #HEC
    ptm <- proc.time()
    out = HEC(X,Y)
    results$HEC$time[j,i] = (proc.time() - ptm)[3]
    if (out >= 0){
      results$HEC$scores[j,i] = 1
    } else{
      results$HEC$scores[j,i] = 0
    }

    #FOM
    ptm <- proc.time()
    out = FOM(X,Y)
    results$FOM$time[j,i] = (proc.time() - ptm)[3]
    if (out >= 0){
      results$FOM$scores[j,i] = 1
    } else{
      results$FOM$scores[j,i] = 0
    }

    # RESIT
    ptm <- proc.time()
    out <- ICML(cbind(normalize01(X),normalize01(Y)), model = train_gp, indtest = indtestHsic, output = FALSE)
    results$ANM$time[j,i] = (proc.time() - ptm)[3]
    if (identical(out$Cd, "->")){
      results$ANM$scores[j,i] = 1
    } else{
      results$ANM$scores[j,i] = 0
    }

    # Direct LiNGAM
    ptm <- proc.time()
    out = LiNGAM_direc(X,Y)
    results$LiNGAM$time[j,i] = (proc.time() - ptm)[3]
    if (out >= 0){
      results$LiNGAM$scores[j,i] = 1
    } else{
      results$LiNGAM$scores[j,i] = 0
    }

    save(file="Res_synth_pairs.RData",results)
    
  }
}


m = rowMeans(results$GRCI$scores[,1:200])
up = m + apply(results$GRCI$scores[,1:200],1,sd)*1.96/sqrt(200)
low = m - apply(results$GRCI$scores[,1:200],1,sd)*1.96/sqrt(200)
print(rbind(m,up,low))


m = rowMeans(results$HEC$scores[,1:200])
up = m + apply(results$HEC$scores[,1:200],1,sd)*1.96/sqrt(200)
low = m - apply(results$HEC$scores[,1:200],1,sd)*1.96/sqrt(200)
print(rbind(m,up,low))

m = rowMeans(results$FOM$scores[,1:200])
up = m + apply(results$FOM$scores[,1:200],1,sd)*1.96/sqrt(200)
low = m - apply(results$FOM$scores[,1:200],1,sd)*1.96/sqrt(200)
print(rbind(m,up,low))

m = rowMeans(results$ANM$scores[,1:200])
up = m + apply(results$ANM$scores[,1:200],1,sd)*1.96/sqrt(200)
low = m - apply(results$ANM$scores[,1:200],1,sd)*1.96/sqrt(200)
print(rbind(m,up,low))

m = rowMeans(results$LiNGAM$scores[,1:200])
up = m + apply(results$LiNGAM$scores[,1:200],1,sd)*1.96/sqrt(200)
low = m - apply(results$LiNGAM$scores[,1:200],1,sd)*1.96/sqrt(200)
print(rbind(m,up,low))


## OVERALL

m = mean(results$GRCI$scores[,1:200])
up = m + sd(results$GRCI$scores[,1:200])*1.96/sqrt(1000)
low = m - sd(results$GRCI$scores[,1:200])*1.96/sqrt(1000)
print(rbind(m,up,low))

m = mean(results$HEC$scores[,1:200])
up = m + sd(results$HEC$scores[,1:200])*1.96/sqrt(1000)
low = m - sd(results$HEC$scores[,1:200])*1.96/sqrt(1000)
print(rbind(m,up,low))

m = mean(results$FOM$scores[,1:200])
up = m + sd(results$FOM$scores[,1:200])*1.96/sqrt(1000)
low = m - sd(results$FOM$scores[,1:200])*1.96/sqrt(1000)
print(rbind(m,up,low))

m = mean(results$ANM$scores[,1:200])
up = m + sd(results$ANM$scores[,1:200])*1.96/sqrt(1000)
low = m - sd(results$ANM$scores[,1:200])*1.96/sqrt(1000)
print(rbind(m,up,low))

m = mean(results$LiNGAM$scores[,1:200])
up = m + sd(results$LiNGAM$scores[,1:200])*1.96/sqrt(1000)
low = m - sd(results$LiNGAM$scores[,1:200])*1.96/sqrt(1000)
print(rbind(m,up,low))

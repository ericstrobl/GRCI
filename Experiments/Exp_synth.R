library(xgboost)
library(pcalg)
library(DirichletReg)
library(Rfast)
library(RANN)
library(treeshap)

ns = c(500,1000,2000)
ps = c(10,30,50)
reps = 500

Gs = lapply(1:reps, function(.) lapply(1:length(ns),function(.) vector("list",length(ps))))
GRCI_res = Gs
ICA_res = Gs
RCI_res = Gs
CO_res = Gs
MS_res = Gs

for (p in 1:length(ps)){
  for (i in 1:reps){
    print(i)
    G = generate_DAG_HNM(p=ps[p],en=2)
    X = sample_DAG_HNM(nsamps = 50000,G)
    X$data[,-G$Y] = normalizeData(X$data[,-G$Y]) # normalize to prevent gaming of variances
    
    #GROUND TRUTH
    Anc = FindAnc(G$graph,G$Y)
    mod = train_xgboost_class(X$E[,Anc,drop=FALSE],X$data[,G$Y])
    unified_mod<- xgboost.unify(mod, X$E[,Anc,drop=FALSE])
    tdata = as.data.frame(X$E[,Anc,drop=FALSE])
    colnames(tdata) = 0:(ncol(tdata)-1)
    
    id = intersect(1:ns[length(ns)],which(X$data[,G$Y]==1))
    shap = matrix(0,ns[length(ns)],ps[p])
    shap[,Anc] <- as.matrix(treeshap(unified_mod, tdata, verbose = FALSE)$shaps[1:ns[length(ns)],])
    shap = shap[,-G$Y]
    
    AncY = Anc - ((1:ncol(G$graph))>G$Y)[Anc] # get ancestor set ready for MSE
    
    # plot(as(G$graph,"graphNEL"))
    # print(G$Y)
    
    for (n in 1:length(ns)){
      idt = intersect(1:ns[n],which(X$data[,G$Y]==1))
      
      Gs[[i]][[n]][[p]] = G
      GT_E = normalizeData_HNM(X$E[1:ns[n],-G$Y])
      
      
      # order and error for GRCI, graph for CO and MS
      ptm <- proc.time()
      outD = learn_DAG_order(X$data[1:ns[n],-G$Y],X$data[1:ns[n],G$Y])
      Gp = outD$G
      time_G = (proc.time() - ptm)[3]
      outL = outD$outL
      GRCI_res[[i]][[n]][[p]]$MSE = eval_L2_HNM(GT_E,outL$X,outL$K,AncY) ## get MSE of errors for GRCI, CO, MS
      
      #GRCI
      ptm <- proc.time()
      out = GRCI(X$data[1:ns[n],-G$Y],X$data[1:ns[n],G$Y],outL)
      time = (proc.time() - ptm)[3]
      GRCI_res[[i]][[n]][[p]]$time = time + outD$time_E
      out$scores = out$scores[idt,]
      GRCI_res[[i]][[n]][[p]]$rank_overlap = eval_scores(shap[idt,],out)
      
      #RCI
      out = RCI(X$data[1:ns[n],-G$Y],X$data[1:ns[n],G$Y])
      RCI_res[[i]][[n]][[p]]$time = out$time_orig
      RCI_res[[i]][[n]][[p]]$time_tree_shap = out$time_tree_shap
      out$scores = out$scores[idt,]
      RCI_res[[i]][[n]][[p]]$rank_overlap = eval_scores(shap[idt,],out)
      out$scores = out$scores_tree_shap[idt,]
      RCI_res[[i]][[n]][[p]]$rank_overlap_tree_shap = eval_scores(shap[idt,],out)
      RCI_res[[i]][[n]][[p]]$MSE = eval_L2_HNM(GT_E,out$E,out$order,AncY)
      
      #ICA
      out = ICA_predict(X$data[1:ns[n],-G$Y],X$data[1:ns[n],G$Y])
      ICA_res[[i]][[n]][[p]]$time = out$time_orig
      ICA_res[[i]][[n]][[p]]$time_tree_shap = out$time_tree_shap
      out$scores = out$scores[idt,]
      ICA_res[[i]][[n]][[p]]$rank_overlap = eval_scores(shap[idt,],out)
      out$scores = out$scores_tree_shap[idt,]
      ICA_res[[i]][[n]][[p]]$rank_overlap_tree_shap = eval_scores(shap[idt,],out)
      ICA_res[[i]][[n]][[p]]$MSE = eval_L2_HNM(GT_E,out$E,out$order,AncY)
      
      #CO
      ptm <- proc.time()
      out = cond_outl_HNM(X$data[1:ns[n],-G$Y],X$data[1:ns[n],G$Y],Gp)
      CO_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3] + time_G
      CO_res[[i]][[n]][[p]]$time_G = time_G
      out$scores = out$scores[idt,]
      CO_res[[i]][[n]][[p]]$rank_overlap = eval_scores(shap[idt,],out)
      
      #MS
      ptm <- proc.time()
      out = model_subst_HNM(X$data[1:ns[n],-G$Y],X$data[1:ns[n],G$Y],Gp,outL$X[,rev(outL$K)])
      MS_res[[i]][[n]][[p]]$time = (proc.time() - ptm)[3] + time_G
      out$scores = matrix(out$scores,length(idt),length(out$scores),byrow=TRUE)
      MS_res[[i]][[n]][[p]]$rank_overlap = eval_scores(shap[idt,],out)
      
      save(file="Results_synth.RData",Gs,GRCI_res,RCI_res,ICA_res,CO_res,MS_res)
      
    }
  }
  
}

## RBO
RBO_GRCI = array(0,c(reps,length(ns),length(ps)))
RBO_RCI = RBO_GRCI
RBO_RCIt = RBO_GRCI
RBO_ICA = RBO_GRCI
RBO_ICAt = RBO_GRCI
RBO_CO = RBO_GRCI
RBO_MS = RBO_GRCI

for (i in 1:reps){
  for (n in 1:length(ns)){
    for (p in 1:length(ps)){
      
      RBO_GRCI[i,n,p] = GRCI_res[[i]][[n]][[p]]$rank_overlap
      RBO_RCI[i,n,p] = RCI_res[[i]][[n]][[p]]$rank_overlap
      RBO_RCIt[i,n,p] = RCI_res[[i]][[n]][[p]]$rank_overlap_tree_shap
      RBO_ICA[i,n,p] = ICA_res[[i]][[n]][[p]]$rank_overlap
      RBO_ICAt[i,n,p] = RCI_res[[i]][[n]][[p]]$rank_overlap_tree_shap
      RBO_CO[i,n,p] = CO_res[[i]][[n]][[p]]$rank_overlap
      RBO_MS[i,n,p] = MS_res[[i]][[n]][[p]]$rank_overlap
      
    }
  }
}

apply(RBO_GRCI,c(2,3),mean)
apply(RBO_RCI,c(2,3),mean)
apply(RBO_RCIt,c(2,3),mean)
apply(RBO_ICA,c(2,3),mean)
apply(RBO_ICAt,c(2,3),mean)
apply(RBO_MS,c(2,3),mean)
apply(RBO_CO,c(2,3),mean)


## MSE
MSE_GRCI = array(0,c(reps,length(ns),length(ps)))
MSE_RCI = MSE_GRCI
MSE_ICA = MSE_GRCI

for (i in 1:reps){
  for (n in 1:3){
    for (p in 1:3){
      
      MSE_GRCI[i,n,p] = GRCI_res[[i]][[n]][[p]]$MSE
      MSE_RCI[i,n,p] = RCI_res[[i]][[n]][[p]]$MSE
      MSE_ICA[i,n,p] = ICA_res[[i]][[n]][[p]]$MSE
      
    }
  }
}

apply(MSE_GRCI,c(2,3),mean)
apply(MSE_RCI,c(2,3),mean)
apply(MSE_ICA,c(2,3),mean)


## time
time_GRCI = array(0,c(reps,length(ns),length(ps)))
time_RCI = time_GRCI
time_RCIt = time_GRCI
time_ICA = time_GRCI
time_ICAt = time_GRCI
time_CO = time_GRCI
time_MS = time_GRCI

for (i in 1:reps){
  for (n in 1:3){
    for (p in 1:3){
      
      time_GRCI[i,n,p] = GRCI_res[[i]][[n]][[p]]$time
      time_RCI[i,n,p] = RCI_res[[i]][[n]][[p]]$time
      time_RCIt[i,n,p] = RCI_res[[i]][[n]][[p]]$time_tree_shap
      time_ICA[i,n,p] = ICA_res[[i]][[n]][[p]]$time
      time_ICAt[i,n,p] = ICA_res[[i]][[n]][[p]]$time_tree_shap
      time_CO[i,n,p] = CO_res[[i]][[n]][[p]]$time
      time_MS[i,n,p] = MS_res[[i]][[n]][[p]]$time
      
    }
  }
}

apply(time_GRCI,c(2,3),mean)
apply(time_RCI,c(2,3),mean)
apply(time_RCIt,c(2,3),mean)
apply(time_ICA,c(2,3),mean)
apply(time_ICAt,c(2,3),mean)
apply(time_MS,c(2,3),mean)
apply(time_CO,c(2,3),mean)


library(xgboost)
library(pcalg)
library(DirichletReg)
library(Rfast)
library(RANN)
library(treeshap)

library(survival)
id = (pbc[,3]!=1) & complete.cases(pbc)
Xo = pbc[id,c(5,11:19)]
X = normalizeData(Xo)
Y = 0+(pbc[id,3]>0) # patient would have died if not transplant, so 0 = alive, 1 = dead or transplant
colnames(X) = c()

reps = 1000

GRCI_res = vector("list",reps)
RCI_res = GRCI_res
ICA_res = GRCI_res
CO_res = GRCI_res
MS_res = GRCI_res

r= length(Y)
for (i in 1:reps){
  print(i)
  
  is = sample(r,r,replace=TRUE)
  Xt = X[is,]
  Yt = Y[is]
  
  idt = which(Yt==1)
  
  ptm <- proc.time()
  outD = learn_DAG_order(Xt,Yt)
  Gp = outD$G
  time_G = (proc.time() - ptm)[3]
  outL = outD$outL
  
  ### RUN ALGORITHMS
  
  #GRCI
  ptm <- proc.time()
  out = GRCI(Xt,Yt,outL)
  time = (proc.time() - ptm)[3]
  GRCI_res[[i]]$time = time + outD$time_E
  out$scores = out$scores[idt,]
  GRCI_res[[i]]$rank_overlap = eval_scores_pbc(out,Xo[is[idt],])
  
  #RCI
  out = RCI(Xt,Yt)
  RCI_res[[i]]$time = out$time_orig
  RCI_res[[i]]$time_tree_shap = out$time_tree_shap
  out$scores = out$scores[idt,]
  RCI_res[[i]]$rank_overlap = eval_scores_pbc(out,Xo[is[idt],])
  out$scores = out$scores_tree_shap[idt,]
  RCI_res[[i]]$rank_overlap_tree_shap = eval_scores_pbc(out,Xo[is[idt],])
  
  out = ICA_predict(Xt,Yt)
  ICA_res[[i]]$time = out$time_orig
  ICA_res[[i]]$time_tree_shap = out$time_tree_shap
  out$scores = out$scores[idt,]
  ICA_res[[i]]$rank_overlap = eval_scores_pbc(out,Xo[is[idt],])
  out$scores = out$scores_tree_shap[idt,]
  ICA_res[[i]]$rank_overlap_tree_shap = eval_scores_pbc(out,Xo[is[idt],])
  
  #CO
  ptm <- proc.time()
  out = cond_outl_HNM(Xt,Yt,Gp)
  CO_res[[i]]$time = (proc.time() - ptm)[3] + time_G
  CO_res[[i]]$time_G = time_G
  out$scores = out$scores[idt,]
  CO_res[[i]]$rank_overlap = eval_scores_pbc(out,Xo[is[idt],])
  
  #MS
  ptm <- proc.time()
  out = model_subst_HNM(Xt,Yt,Gp,outL$X[,rev(outL$K)])
  MS_res[[i]]$time = (proc.time() - ptm)[3] + time_G
  out$scores = matrix(out$scores,length(is),length(out$scores),byrow=TRUE)
  out$scores = out$scores[idt,]
  MS_res[[i]]$rank_overlap = eval_scores_pbc(out,Xo[is[idt],])
  
  
  ### SAVE RESULTS
  save(file="Res_PBC2_time.RData",GRCI_res,RCI_res,ICA_res,CO_res,MS_res)
  
}


### PBC RESULTS
RBO = matrix(0,reps,7)
for (i in 1:reps){
  RBO[i,1] = GRCI_res[[i]]$rank_overlap
  RBO[i,2] = RCI_res[[i]]$rank_overlap
  RBO[i,3] = RCI_res[[i]]$rank_overlap_tree_shap
  RBO[i,4] = ICA_res[[i]]$rank_overlap
  RBO[i,5] = ICA_res[[i]]$rank_overlap_tree_shap
  RBO[i,6] = CO_res[[i]]$rank_overlap
  RBO[i,7] = MS_res[[i]]$rank_overlap
}

print(colMeans(RBO))
print(apply(RBO,2,sd)/sqrt(reps)*1.96)


### PBC TIMING
## time
times = matrix(0,reps,7)
for (i in 1:reps){
  times[i,1] = GRCI_res[[i]]$time
  times[i,2] = RCI_res[[i]]$time
  times[i,3] = RCI_res[[i]]$time_tree_shap
  times[i,4] = ICA_res[[i]]$time
  times[i,5] = ICA_res[[i]]$time_tree_shap
  times[i,6] = CO_res[[i]]$time
  times[i,7] = MS_res[[i]]$time
}

print(colMeans(times))
print(apply(times,2,sd)/sqrt(reps)*1.96)

ICA_predict <- function(X,Y){
  
  time = proc.time()
  require(ica)
  c = ncol(X)
  time = proc.time()
  mm = icafast(X, c)
  
  require(RcppHungarian)
  
  id = HungarianSolver(t(1/abs(mm$M)))$pairs[,2]
  E = mm$S[,order(id)]
  
  time_search = (proc.time() - time)[3]
  
  # mm$W = t(mm$W)
  # id = c()
  # for (i in 1:c){
  #   Mi = abs(mm$W[,i]) # what is the variable in X that contributes most to source S[,i]
  #   id  = c(id,which(Mi == max(Mi))[1])
  # }
  # E = mm$S[,order(id)]
  
  time = proc.time()
  require(randomForest)
  
  imp = randomForest(E, y=factor(Y), localImp=TRUE)$localImp
  # imp = randomForest(E, y=factor(Y), importance=TRUE)$importance
  # print(imp)
  time_RF = (proc.time() - time)[3]
  
  
  time = proc.time()
  require(xgboost)
  mod = train_xgboost_class(E,Y)
  
  ## COMPUTE SHAPLEY VALUES
  # print('Shapley values')
  require(treeshap)
  unified_mod <- xgboost.unify(mod, E)
  tdata = as.data.frame(E)
  colnames(tdata) = 0:(ncol(tdata)-1)
  shap <- as.matrix(treeshap(unified_mod, tdata, verbose = FALSE)$shaps)
  time_tree_shap = (proc.time() - time)[3]
  
  return(list(scores = t(imp), scores_tree_shap = shap, order = 1:c, E = E, 
              time_orig = time_search + time_RF, time_tree_shap = time_search + time_tree_shap))
  
}



train_xgboost_class <- function(X,y,pos=FALSE){
  n <- length(y)
  pars=list()
  pars$nrounds <- 50
  pars$max_depth <- c(1,3,4,5,6)
  pars$CV.folds <- 10
  pars$ncores <- 1
  pars$early_stopping <- 10
  pars$silent <- TRUE
  
  
  X <- as.matrix(X)
  
  num.folds <- pars$CV.folds
  rmse <- matrix(0,pars$nrounds, length(pars$max_depth))
  whichfold <- sample(rep(1:num.folds, length.out = n))
  for(j in 1:length(pars$max_depth)){
    max_depth <- pars$max_depth[j]
    for(i in 1:num.folds){
      dtrain <- xgb.DMatrix(data = data.matrix(X[whichfold != i,]), label=y[whichfold != i])
      dtest <- xgb.DMatrix(data = data.matrix(X[whichfold == i,]), label=y[whichfold == i])
      watchlist <- list(train = dtrain, test = dtest)
      
      if(pars$ncores > 1){
        bst <- xgb.train(data = dtrain, objective="binary:logistic", eval_metric = "logloss", nthread = pars$ncores, watchlist = watchlist, nrounds = pars$nrounds, max_depth = max_depth, verbose = FALSE, early_stopping_rounds = pars$early_stopping, callbacks = list(cb.evaluation.log()))
      } else{
        bst <- xgb.train(data = dtrain, objective="binary:logistic", eval_metric = "logloss", nthread = 1, watchlist = watchlist, nrounds = pars$nrounds, max_depth = max_depth, verbose = FALSE, early_stopping_rounds = pars$early_stopping, callbacks = list(cb.evaluation.log()))
      }
      # print(bst$evaluation_log)
      newscore <- (bst$evaluation_log$test_logloss)^1
      # print(bst$evaluation_log$test_rmse)
      if(length(newscore) < pars$nrounds){
        newscore <- c(newscore, rep(Inf, pars$nrounds - length(newscore)))
      }
      rmse[,j] <- rmse[,j] + newscore
    }
  }
  
  mins <- arrayInd(which.min(rmse),.dim = dim(rmse))
  final.nrounds <- mins[1]
  final.max_depth <- pars$max_depth[mins[2]]
  
  
  dtrain <- xgb.DMatrix(data = data.matrix(X), label=y)
  bstY <- xgb.train(data = dtrain, nrounds = final.nrounds, max_depth = final.max_depth,
                    verbose = !pars$silent)
  
  return(bstY)
}
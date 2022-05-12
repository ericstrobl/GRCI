# GRCI <- function(X,Y){
GRCI <- function(X,Y,outL=NULL){
  
  if (is.null(outL)){
    X = normalizeData(X)

    ## SKELETON DISCOVERY
    # print('skeleton')
    suffStat = list()
    suffStat$data = X
    G = pc(suffStat, earth_wrap, alpha=0.10, p=ncol(X))
    G = as(G@graph, "matrix")
    G = ((G + t(G))>0)
    # print(G)

    ## EXTRACT ERRORS
    # print('errors')
    outL = DirectHNM_fast_Y(X,Y,G) # Local Plus
  }
  
  ## LOGISTIC REGRESSION
  # print('predictive model')
  require(xgboost)
  mod = train_xgboost_class(outL$X,Y)
  
  ## COMPUTE SHAPLEY VALUES
  # print('Shapley values')
  require(treeshap)
  unified_mod <- xgboost.unify(mod, outL$X)
  tdata = as.data.frame(outL$X)
  colnames(tdata) = 0:(ncol(tdata)-1)
  shap <- as.matrix(treeshap(unified_mod, tdata, verbose = FALSE)$shaps)
  
  return(list(E=outL$X, order=outL$K, scores=shap, G=G))
}

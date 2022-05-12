#' Root Causal Inference (RCI) algorithm.
#'
#' @param X a matrix of continuous variables, rows are samples and columns are variables
#' @param Y a vector of a binary target
#' @param alpha the alpha value for the t-test, default is 0.2
#' @return A list containing Shapley values \code{scores}, the errors \code{E}, the variable ordering \code{K}, and time of Local Plus \code{time}
#' @export

RCI <- function(X,Y,alpha=0.2){
  
  time = proc.time()
  outL = DirectLiNGAM_fast_Y(X,Y) # Local Plus
  # outL = DirectLiNGAM(X)
  K = outL$K;
  if (length(K)==0){
    delta = matrix(0,nrow(X),ncol(X))
    return(list(delta=delta,K=K))
  }
  E = outL$X[,order(K),drop=FALSE]
  K = sort(K)
  time_search = (proc.time() - time)[3]
  
  
  time = proc.time()
  ## patient-specific root causes
  beta = glm.fit(cbind(E,1),Y,family=binomial())$coefficients[1:length(K)]
  # print(beta)
  if (length(beta)==1){
    scores = E * beta
  } else{
    scores = E %*% diag(beta)
  }
  time_linear = (proc.time() - time)[3]
  
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
  
  colnames(scores) = K
  colnames(E) = K
  
  return(list(scores=scores, order=K, E=E, scores_tree_shap = shap,
              time_orig = time_search + time_linear, 
              time_tree_shap = time_search + time_tree_shap))
}

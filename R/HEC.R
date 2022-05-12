polyfit <- function(x,y,n_glob){ ###
  n_loc = length(x)
  opt_cost = 0
  res = 0
  for (i in 1:3){
    m = lm.fit(cbind(1,polym(x,degree=i)),y)
    # coeff = m$coefficients
    res = sum(m$residuals^2)
    
    res = res + 0.01
    data_cost = log2(res/n_loc) * n_loc
    model_cost = (i+2)*log2(n_glob)
    cost = data_cost + model_cost
    if ( (i == 1) | (cost < opt_cost) ){
      opt_cost = cost
      # opt_coeff = coeff
      # opt_deg = i
    }
  }
  
  return(opt_cost)
}

precomputeScores <- function(binsx, binsy, fitting_fun){ 
  beta = length(binsx)
  scores = matrix(0,beta,beta)
  
  # model_params = lapply(1:beta, function(.) vector("list", beta))
  n = length(unlist(binsx))
  for (i in 0:(beta-1)){
      start = 0
      stop = start + i + 1
      while (stop <= beta){
        data_split_x = unlist(binsx[(start+1):stop])
        data_split_y = unlist(binsy[(start+1):stop])
        score_fit = fitting_fun(data_split_x, data_split_y, n)
        # deg = fit$opt_deg
        # coeff = fit$opt_coeff
        
        # model_params[[i]][[start]]$deg = deg
        # model_params[[i]][[start]]$coeff = coeff
        
        scores[i+1, start+1] = score_fit
        start = start + 1
        stop = stop + 1
      }
  }
  
  return( list(scores=scores, beta=beta) )
  
}

HEC_Opt <- function(binned_x, binned_y, fitting_fun=polyfit, show=FALSE){
  pS = precomputeScores(binned_x, binned_y, fitting_fun)
  precomputed_scores=pS$scores
  num_bins = pS$beta
  
  opt_cost = precomputed_scores
  
  # opt_split = array(0, dim = c(num_bins, num_bins, 4))
  # opt_split[1, ,] = -1
  
  for (i in seq_len(num_bins-1)){
    start = 0
    stop = start + i + 1
    while (stop <= num_bins){
      
      min_cost = opt_cost[i+1, start+1]
      # opt_split[i, start,] = -1
      for (j in 0:(i-1)){
        split1cost = opt_cost[j+1, start+1]
        split2size = i - j - 1
        split2loc = start + j + 1
        split2cost = opt_cost[split2size+1, split2loc+1]
        new_cost = split1cost + split2cost
        
        if (min_cost > new_cost){
          min_cost = new_cost
          # opt_split[i, start, ] = c(j, start, split2size, split2loc)
        }

      }
      opt_cost[i+1, start+1] = min_cost
      start = start + 1
      stop = stop + 1
    }
  }
  
  final_cost = opt_cost[num_bins, 1]
  return(final_cost)
  
  
}

HEC <- function(x,y){

  list1 = widthBinning(x,y,beta=20)
  x1 = list1$binsx
  y1 = list1$binsy
  
  list2 = widthBinning(y,x,beta=20)
  x2 = list2$binsx
  y2 = list2$binsy

  score1 = HEC_Opt(x1, y1, polyfit)
  score2 = HEC_Opt(x2, y2, polyfit)
  
  if (score1 <= score2){
    return(1)
  } else{
    return(-1)
  }

}

widthBinning <- function(x, y, beta=20, min_support=10){ ###
  ind = order(x)
  y = y[ind]
  x = x[ind]
  
  y = normalize01(y)
  x = normalize01(x)
  
  step = (max(x) - min(x)) / beta  # 0.05 step described in paper
  threshold = min(x) + step
  count = 1
  
  bin_indices = c()

  for (index in 0:(length(x)-1)){
    if (x[index+1] <= threshold){
      next
    } else{
      if (length(bin_indices) == 0){
        last = 0
      } else{
        last = bin_indices[length(bin_indices)]
      }
      check = x[(last+1):index]
      support = length(unique(check))
      
      if (support >= min_support){
        bin_indices = c(bin_indices,index)
      }
      
      count = count + 1
      threshold = threshold + step
      
      if (count == beta){
        break
      }
    }
  }
  
  rest = x[(bin_indices[length(bin_indices)]+1):length(x)]
  # support = length(unique(rest))
  # 
  # if (support < min_support){
  #   bin_indices = bin_indices[-length(bin_indices)]
  # }
  
  bin_indices = c(bin_indices[-length(bin_indices)],bin_indices[length(bin_indices)]+length(rest))
  
  binsx = list()
  s = 0
  for (b in 1:length(bin_indices)){
    binsx[[b]] = x[(s+1):bin_indices[b]]
    s = bin_indices[b]
  }
  
  binsy = list()
  s = 0
  for (b in 1:length(bin_indices)){
    binsy[[b]] = y[(s+1):bin_indices[b]]
    s = bin_indices[b]
  }
  
  return(list( binsx = binsx, binsy = binsy, beta = beta ) )
  
}

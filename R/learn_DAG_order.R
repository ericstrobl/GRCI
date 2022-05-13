learn_DAG_order <- function(X,Y){
  
  ptm <- proc.time()
  X = normalizeData_HNM(X)
  
  suffStat = list()
  suffStat$data = X
  p = ncol(X)
  G = my_skeleton(matrix(TRUE,p,p),suffStat, earth_wrap, alpha=0.10, p=p)
  G = as(G@graph, "matrix")
  G = ((G + t(G))>0)
  # print(G)
  
  ## extract errors and order
  outL = DirectHNM_fast_Y(X,Y,G) # Local Plus, extracts reverse partial order
  time_E = (proc.time() - ptm)[3]
  
  ## use order to find skeleton and direct edges
  outL$K = rev(outL$K) # partial order
  suffStat$data = cbind(X, as.matrix(Y))
  G[outL$K,outL$K] = upper.tri(matrix(0,p,p),diag=FALSE)
  G = rbind(cbind(G,TRUE),FALSE) # because Y is always last row and last column
  
  G = my_skeleton(G,suffStat, earth_wrap, alpha=0.05, p=(p+1)) # find skeleton
  G = as(G@graph, "matrix")
  Gt = G[outL$K,outL$K]
  Gt[!upper.tri(Gt)]=FALSE
  G[outL$K,outL$K] = Gt[1:p,1:p]
  G[nrow(G),]=FALSE
  
  outL$K = rev(outL$K) #revert back to reverse partial order
  return(list(G=G,outL=outL,time_E=time_E))
}

my_skeleton <- function (G_old, suffStat, indepTest, alpha, labels, p, method = "stable",
                         m.max = Inf, fixedGaps = NULL, 
                         fixedEdges = NULL, NAdelete = TRUE, numCores = 1, verbose = FALSE) 
{
  cl <- match.call()
  if (!missing(p)) 
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) 
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) 
      p <- length(labels)
    else if (p != length(labels)) 
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedGaps), c(p, p))) 
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps))) 
    stop("fixedGaps must be symmetric")
  else G <- !fixedGaps
  diag(G) <- FALSE
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p))) 
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges))) 
    stop("fixedEdges must be symmetric")
  stopifnot((is.integer(numCores) || is.numeric(numCores)) && 
              numCores > 0)
  
  pval <- NULL
  sepset <- lapply(seq_p, function(.) vector("list", 
                                             p))
  pMax <- matrix(-Inf, nrow = p, ncol = p)
  diag(pMax) <- 1
  done <- FALSE
  ord <- 0L
  n.edgetests <- numeric(1)
  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord + 1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE) ####
    ind <- ind[order(ind[, 1]), ]
    remEdges <- nrow(ind)
    
    G.l <- split(G, gl(p, p))
    
    for (i in 1:remEdges) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      if (G[y, x] && !fixedEdges[y, x]) {
        nbrsBool <- G.l[[x]] ####
        nbrsBool[y] <- FALSE
        nbrs <- seq_p[nbrsBool]
        nbrs <- intersect(nbrs, which(G_old[,x]>0)) ###
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) 
            done <- FALSE
          S <- seq_len(ord)
          repeat {
            n.edgetests[ord1] <- n.edgetests[ord1] + 
              1
            pval <- indepTest(x, y, nbrs[S], suffStat)
            if (is.na(pval)) 
              pval <- as.numeric(NAdelete)
            if (pMax[x, y] < pval) 
              pMax[x, y] <- pval
            if (pval >= alpha) {
              G[x, y] <- G[y, x] <- FALSE
              sepset[[x]][[y]] <- nbrs[S]
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, 
                                    S)
              if (nextSet$wasLast) 
                break
              S <- nextSet$nextSet
            }
          }
        }
      }
    }
    ord <- ord + 1L
  }
  for (i in 1:(p - 1)) {
    for (j in 2:p) pMax[i, j] <- pMax[j, i] <- max(pMax[i, 
                                                        j], pMax[j, i])
  }
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  }
  else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL")
  }
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}

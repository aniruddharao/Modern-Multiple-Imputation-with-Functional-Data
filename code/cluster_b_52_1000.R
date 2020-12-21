library(ggplot2)
library(reshape2)
library(MASS)
library(mice)
library(lattice)
library(missForest)
library(fdapace)
library(refund)
library(Metrics)
library(grf)



pacemissForest <- function(xmis, maxiter = 10, ntree = 100, variablewise = FALSE,
                           decreasing = FALSE, verbose = FALSE,
                           mtry = floor(sqrt(ncol(xmis))), replace = TRUE,
                           classwt = NULL, cutoff = NULL, strata = NULL,
                           sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                           xtrue = NA, parallelize = c('no', 'variables', 'forests'))
{ 
  n <- nrow(xmis)
  p <- ncol(xmis)
  if (!is.null(classwt))
    stopifnot(length(classwt) == p, typeof(classwt) == 'list')
  if (!is.null(cutoff))
    stopifnot(length(cutoff) == p, typeof(cutoff) == 'list')
  if (!is.null(strata))
    stopifnot(length(strata) == p, typeof(strata) == 'list')
  if (!is.null(nodesize))
    stopifnot(length(nodesize) == 2)
  
  ## remove completely missing variables
  if (any(apply(is.na(xmis), 2, sum) == n)){
    indCmis <- which(apply(is.na(xmis), 2, sum) == n)
    xmis <- xmis[,-indCmis]
    p <- ncol(xmis)
    cat('  removed variable(s)', indCmis,
        'due to the missingness of all entries\n')
  } 
  
  ## return feedback on parallelization setup
  parallelize <- match.arg(parallelize)
  if (parallelize %in% c('variables', 'forests')) {
    if (getDoParWorkers() == 1) {
      stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
    } else if (verbose) {
      if (parallelize == 'variables') {
        cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
      } else {
        cat("  parallelizing computation of the random forest model objects\n")
      }
    }
    if (getDoParWorkers() > p){
      stop('The number of parallel cores should not exceed the number of variables (p=', p, ")")
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  ## perform initial S.W.A.G. on xmis (mean imputation)
  
  ximp <- xmis
  xAttrib <- lapply(xmis, attributes)
  
  p=m
  res= FPCA(mlisty, mlistx,
            list(dataType='Sparse', error=FALSE,nRegGrid = m, kernel='epan', verbose=TRUE))
  
  X.pace=fitted.values(res)
  varType <- character(p)
  for (t.co in 1:p){
    if (is.null(xAttrib[[t.co]])){
      varType[t.co] <- 'numeric'
      ximp[is.na(xmis[,t.co]),t.co] <- X.pace[is.na(xmis[,t.co]),t.co]
    } else {
      varType[t.co] <- 'factor'
      ## take the level which is more 'likely' (majority vote)
      max.level <- max(table(ximp[,t.co]))
      ## if there are several classes which are major, sample one at random
      class.assign <- sample(names(which(max.level == summary(ximp[,t.co]))), 1)
      ## it shouldn't be the NA class
      if (class.assign != "NA's"){
        ximp[is.na(xmis[,t.co]),t.co] <- class.assign
      } else {
        while (class.assign == "NA's"){
          class.assign <- sample(names(which(max.level ==
                                               summary(ximp[,t.co]))), 1)
        }
        ximp[is.na(xmis[,t.co]),t.co] <- class.assign
      }
    }
  }
  
  ## extract missingness pattern
  NAloc <- is.na(xmis)            # where are missings
  noNAvar <- apply(NAloc, 2, sum) # how many are missing in the vars
  sort.j <- order(noNAvar)        # indices of increasing amount of NA in vars
  if (decreasing)
    sort.j <- rev(sort.j)
  sort.noNAvar <- noNAvar[sort.j]
  
  ## compute a list of column indices for variable parallelization
  nzsort.j <- sort.j[sort.noNAvar > 0]
  if (parallelize == 'variables') {
    '%cols%' <- get('%dopar%')
    idxList <- as.list(isplitVector(nzsort.j, chunkSize=getDoParWorkers()))
  } 
  #   else {
  #     ## force column loop to be sequential
  #     '%cols%' <- get('%do%')
  #     idxList <- nzsort.j
  #   }
  
  ## output
  Ximp <- vector('list', maxiter)
  
  ## initialize parameters of interest
  iter <- 0
  k <- length(unique(varType))
  convNew <- rep(0, k)
  convOld <- rep(Inf, k)
  OOBerror <- numeric(p)
  names(OOBerror) <- varType
  
  ## setup convergence variables w.r.t. variable types
  if (k == 1){
    if (unique(varType) == 'numeric'){
      names(convNew) <- c('numeric')
    } else {
      names(convNew) <- c('factor')
    }
    convergence <- c()
    OOBerr <- numeric(1)
  } else {
    names(convNew) <- c('numeric', 'factor')
    convergence <- matrix(NA, ncol = 2)
    OOBerr <- numeric(2)
  }
  
  ## function to yield the stopping criterion in the following 'while' loop
  stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
    k <- length(unique(varType))
    if (k == 1){
      (convNew < convOld) & (iter < maxiter)
    } else {
      ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)
    }
  }
  
  ## iterate missForest
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)){
    if (iter != 0){
      convOld <- convNew
      OOBerrOld <- OOBerr
    }
    cat("  missForest iteration", iter+1, "in progress...")
    t.start <- proc.time()
    ximp.old <- ximp
    
    if (parallelize=="variables"){
      for (idx in idxList) {
        results <- foreach(varInd=idx, .packages='randomForest') %cols% {
          obsi <- !NAloc[,varInd] # which i's are observed
          misi <- NAloc[,varInd] # which i's are missing
          obsY <- ximp[obsi, varInd] # training response
          obsX <- ximp[obsi, seq(1, p)[-varInd]] # training variables
          misX <- ximp[misi, seq(1, p)[-varInd]] # prediction variables
          typeY <- varType[varInd]
          if (typeY == 'numeric'){
            RF <- randomForest(
              x = obsX,
              y = obsY,
              ntree = ntree,
              mtry = mtry,
              replace = replace,
              sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
              nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
              maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
            ## record out-of-bag error
            oerr <- RF$mse[ntree]
            #           }
            ## predict missing values in column varInd
            misY <- predict(RF, misX)
          } else { # if Y is categorical          
            obsY <- factor(obsY) ## remove empty classes
            summarY <- summary(obsY)
            if (length(summarY) == 1){ ## if there is only one level left
              oerr <- 0
              misY <- factor(rep(names(summarY), length(misi)))
            } else {
              RF <- randomForest(
                x = obsX,
                y = obsY,
                ntree = ntree,
                mtry = mtry,
                replace = replace,
                classwt = if (!is.null(classwt)) classwt[[varInd]] else
                  rep(1, nlevels(obsY)),
                cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else
                  rep(1/nlevels(obsY), nlevels(obsY)),
                strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                  if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
              ## record out-of-bag error
              oerr <- RF$err.rate[[ntree,1]]
              #             }
              ## predict missing values in column varInd
              misY <- predict(RF, misX)
            }
          }
          list(varInd=varInd, misY=misY, oerr=oerr)
        }
        ## update the master copy of the data
        for (res in results) {
          misi <- NAloc[,res$varInd]
          ximp[misi, res$varInd] <- res$misY
          OOBerror[res$varInd] <- res$oerr
        }
      }
    } else { # if parallelize != "variables"
      for (s in 1:p) {
        varInd <- sort.j[s]
        if (noNAvar[[varInd]] != 0) {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- ximp[obsi, varInd]
          obsX <- ximp[obsi, seq(1, p)[-varInd]]
          misX <- ximp[misi, seq(1, p)[-varInd]]
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            if (parallelize == 'forests') {
              xntree <- NULL
              RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                            .combine='combine', .multicombine=TRUE,
                            .packages='randomForest') %dopar% {
                              randomForest( x = obsX,
                                            y = obsY,
                                            ntree = xntree,
                                            mtry = mtry,
                                            replace = replace,
                                            sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                              if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                                            nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                                            maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                            }
              ## record out-of-bag error
              OOBerror[varInd] <- mean((predict(RF) - RF$y) ^ 2, na.rm=TRUE)
              #               OOBerror[varInd] <- RF$mse[ntree]
            } else {
              RF <- randomForest( x = obsX,
                                  y = obsY,
                                  ntree = ntree,
                                  mtry = mtry,
                                  replace = replace,
                                  sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                    if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                                  nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
              ## record out-of-bag error
              OOBerror[varInd] <- RF$mse[ntree]
            }
            misY <- predict(RF, misX)
          } else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              misY <- factor(rep(names(summarY), sum(misi)))
            } else {
              if (parallelize == 'forests') {
                RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                              .combine='combine', .multicombine=TRUE,
                              .packages='randomForest') %dopar% {
                                randomForest(
                                  x = obsX,
                                  y = obsY,
                                  ntree = xntree,
                                  mtry = mtry,
                                  replace = replace,
                                  classwt = if (!is.null(classwt)) classwt[[varInd]] else
                                    rep(1, nlevels(obsY)),
                                  cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else
                                    rep(1/nlevels(obsY), nlevels(obsY)),
                                  strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                                  sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                    if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                                  nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                              }
                ## record out-of-bag error
                ne <- as.integer(predict(RF)) != as.integer(RF$y)
                ne <- ne[! is.na(ne)]
                OOBerror[varInd] <- sum(ne) / length(ne)
              } else {
                RF <- randomForest(x = obsX, 
                                   y = obsY, 
                                   ntree = ntree, 
                                   mtry = mtry, 
                                   replace = replace, 
                                   classwt = if (!is.null(classwt)) classwt[[varInd]] else 
                                     rep(1, nlevels(obsY)),
                                   cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else 
                                     rep(1/nlevels(obsY), nlevels(obsY)),
                                   strata = if (!is.null(strata)) strata[[varInd]] else obsY, 
                                   sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else 
                                     if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)), 
                                   nodesize = if (!is.null(nodesize)) nodesize[2] else 5, 
                                   maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                ## record out-of-bag error
                OOBerror[varInd] <- RF$err.rate[[ntree, 1]]
              }
              ## predict missing parts of Y
              misY <- predict(RF, misX)
            }
          }
          ximp[misi, varInd] <- misY
        }
      }
    }
    cat('done!\n')
    
    iter <- iter+1
    Ximp[[iter]] <- ximp
    
    t.co2 <- 1
    ## check the difference between iteration steps
    for (t.type in names(convNew)){
      t.ind <- which(varType == t.type)
      if (t.type == "numeric"){
        convNew[t.co2] <- sum((ximp[,t.ind]-ximp.old[,t.ind])^2)/sum(ximp[,t.ind]^2)
      } else {
        dist <- sum(as.character(as.matrix(ximp[,t.ind])) != as.character(as.matrix(ximp.old[,t.ind])))
        convNew[t.co2] <- dist / (n * sum(varType == 'factor'))
      }
      t.co2 <- t.co2 + 1
    }
    
    ## compute estimated imputation error
    if (!variablewise){
      NRMSE <- sqrt(mean(OOBerror[varType=='numeric'])/
                      var(as.vector(as.matrix(xmis[,varType=='numeric'])),
                          na.rm = TRUE))
      PFC <- mean(OOBerror[varType=='factor'])
      if (k==1){
        if (unique(varType)=='numeric'){
          OOBerr <- NRMSE
          names(OOBerr) <- 'NRMSE'
        } else {
          OOBerr <- PFC
          names(OOBerr) <- 'PFC'
        }
      } else {
        OOBerr <- c(NRMSE, PFC)
        names(OOBerr) <- c('NRMSE', 'PFC')
      }
    } else {
      OOBerr <- OOBerror
      names(OOBerr)[varType=='numeric'] <- 'MSE'
      names(OOBerr)[varType=='factor'] <- 'PFC'
    }
    
    if (any(!is.na(xtrue))){
      err <- suppressWarnings(mixError(ximp, xmis, xtrue))
    }
    
    ## return status output, if desired
    if (verbose){
      delta.start <- proc.time() - t.start
      if (any(!is.na(xtrue))){
        cat("    error(s):", err, "\n")
      }
      cat("    estimated error(s):", OOBerr, "\n")
      cat("    difference(s):", convNew, "\n")
      cat("    time:", delta.start[3], "seconds\n\n")
    }
  }#end while((convNew<convOld)&(iter<maxiter)){
  
  ## produce output w.r.t. stopping rule
  if (iter == maxiter){
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr)
    } else {
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, error = err)
    }
  } else {
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld)
    } else {
      out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld,
                  error = suppressWarnings(mixError(Ximp[[iter-1]], xmis, xtrue)))
    }
  }
  class(out) <- 'missForest'
  return(out)
}




missForest1 <- function(xmis, maxiter = 10, ntree = 100, variablewise = FALSE,
                        decreasing = FALSE, verbose = FALSE,
                        mtry = floor(sqrt(ncol(xmis))), replace = TRUE,
                        classwt = NULL, cutoff = NULL, strata = NULL,
                        sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                        xtrue = NA, parallelize = c('no', 'variables', 'forests'))
{ ## ----------------------------------------------------------------------
  ## Arguments:
  ## xmis         = data matrix with missing values
  ## maxiter      = stop after how many iterations (default = 10)
  ## ntree        = how many trees are grown in the forest (default = 100)
  ## variablewise = (boolean) return OOB errors for each variable separately
  ## decreasing   = (boolean) if TRUE the columns are sorted with decreasing
  ##                amount of missing values
  ## verbose      = (boolean) if TRUE then missForest returns error estimates,
  ##                runtime and if available true error during iterations
  ## mtry         = how many variables should be tried randomly at each node
  ## replace      = (boolean) if TRUE bootstrap sampling (with replacements)
  ##                is performed, else subsampling (without replacements)
  ## classwt      = list of priors of the classes in the categorical variables
  ## cutoff       = list of class cutoffs for each categorical variable
  ## strata       = list of (factor) variables used for stratified sampling
  ## sampsize     = list of size(s) of sample to draw
  ## nodesize     = minimum size of terminal nodes, vector of length 2, with
  ##                number for continuous variables in the first entry and
  ##                number for categorical variables in the second entry
  ## maxnodes     = maximum number of terminal nodes for individual trees
  ## xtrue        = complete data matrix
  ##
  ## ----------------------------------------------------------------------
  ## Author: Daniel Stekhoven, stekhoven@stat.math.ethz.ch
  
  ## stop in case of wrong inputs passed to randomForest
  n <- nrow(xmis)
  
  p <- ncol(xmis)
  if (!is.null(classwt))
    stopifnot(length(classwt) == p, typeof(classwt) == 'list')
  if (!is.null(cutoff))
    stopifnot(length(cutoff) == p, typeof(cutoff) == 'list')
  if (!is.null(strata))
    stopifnot(length(strata) == p, typeof(strata) == 'list')
  if (!is.null(nodesize))
    stopifnot(length(nodesize) == 2)
  
  ## remove completely missing variables
  if (any(apply(is.na(xmis), 2, sum) == n)){
    indCmis <- which(apply(is.na(xmis), 2, sum) == n)
    xmis <- xmis[,-indCmis]
    p <- ncol(xmis)
    cat('  removed variable(s)', indCmis,
        'due to the missingness of all entries\n')
  } 
  
  ## return feedback on parallelization setup
  parallelize <- match.arg(parallelize)
  if (parallelize %in% c('variables', 'forests')) {
    if (getDoParWorkers() == 1) {
      stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
    } else if (verbose) {
      if (parallelize == 'variables') {
        cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
      } else {
        cat("  parallelizing computation of the random forest model objects\n")
      }
    }
    if (getDoParWorkers() > p){
      stop('The number of parallel cores should not exceed the number of variables (p=', p, ")")
    }
  }
  
  ## perform initial S.W.A.G. on xmis (mean imputation)
  ximp <- xmis
  xAttrib <- lapply(xmis, attributes)
  varType <- character(p)
  for (t.co in 1:p){
    if (is.null(xAttrib[[t.co]])){
      varType[t.co] <- 'numeric'
      ximp[is.na(xmis[,t.co]),t.co] <- mean(xmis[,t.co], na.rm = TRUE)
    } else {
      varType[t.co] <- 'factor'
      ## take the level which is more 'likely' (majority vote)
      max.level <- max(table(ximp[,t.co]))
      ## if there are several classes which are major, sample one at random
      class.assign <- sample(names(which(max.level == summary(ximp[,t.co]))), 1)
      ## it shouldn't be the NA class
      if (class.assign != "NA's"){
        ximp[is.na(xmis[,t.co]),t.co] <- class.assign
      } else {
        while (class.assign == "NA's"){
          class.assign <- sample(names(which(max.level ==
                                               summary(ximp[,t.co]))), 1)
        }
        ximp[is.na(xmis[,t.co]),t.co] <- class.assign
      }
    }
  }
  
  ## extract missingness pattern
  NAloc <- is.na(xmis)            # where are missings
  noNAvar <- apply(NAloc, 2, sum) # how many are missing in the vars
  sort.j <- order(noNAvar)        # indices of increasing amount of NA in vars
  if (decreasing)
    sort.j <- rev(sort.j)
  sort.noNAvar <- noNAvar[sort.j]
  
  ## compute a list of column indices for variable parallelization
  nzsort.j <- sort.j[sort.noNAvar > 0]
  if (parallelize == 'variables') {
    '%cols%' <- get('%dopar%')
    idxList <- as.list(isplitVector(nzsort.j, chunkSize=getDoParWorkers()))
  } 
  #   else {
  #     ## force column loop to be sequential
  #     '%cols%' <- get('%do%')
  #     idxList <- nzsort.j
  #   }
  
  ## output
  Ximp <- vector('list', maxiter)
  
  ## initialize parameters of interest
  iter <- 0
  k <- length(unique(varType))
  convNew <- rep(0, k)
  convOld <- rep(Inf, k)
  OOBerror <- numeric(p)
  names(OOBerror) <- varType
  
  ## setup convergence variables w.r.t. variable types
  if (k == 1){
    if (unique(varType) == 'numeric'){
      names(convNew) <- c('numeric')
    } else {
      names(convNew) <- c('factor')
    }
    convergence <- c()
    OOBerr <- numeric(1)
  } else {
    names(convNew) <- c('numeric', 'factor')
    convergence <- matrix(NA, ncol = 2)
    OOBerr <- numeric(2)
  }
  
  ## function to yield the stopping criterion in the following 'while' loop
  stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
    k <- length(unique(varType))
    if (k == 1){
      (convNew < convOld) & (iter < maxiter)
    } else {
      ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)
    }
  }
  
  ## iterate missForest
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)){
    if (iter != 0){
      convOld <- convNew
      OOBerrOld <- OOBerr
    }
    cat("  missForest iteration", iter+1, "in progress...")
    t.start <- proc.time()
    ximp.old <- ximp
    
    if (parallelize=="variables"){
      for (idx in idxList) {
        results <- foreach(varInd=idx, .packages='grf') %cols% {
          obsi <- !NAloc[,varInd] # which i's are observed
          misi <- NAloc[,varInd] # which i's are missing
          obsY <- ximp[obsi, varInd] # training response
          obsX <- ximp[obsi, seq(1, p)[-varInd]] # training variables
          misX <- ximp[misi, seq(1, p)[-varInd]] # prediction variables
          typeY <- varType[varInd]
          if (typeY == 'numeric'){
            RF <- local_linear_forest(
              X = obsX,
              Y = obsY,
              sample.fraction = 0.5, mtry = NULL,
              num.trees = 2000, num.threads = NULL, min.node.size = NULL,
              honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
              alpha = NULL, imbalance.penalty = NULL,
              compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
              samples_per_cluster = NULL, tune.parameters = F,
              num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
            ## record out-of-bag error
            oerr <- 0 #RF$mse[num.trees]
            #           }
            ## predict missing values in column varInd
            misY <- predict(RF, misX)$predictions
          } else { # if Y is categorical          
            obsY <- factor(obsY) ## remove empty classes
            summarY <- summary(obsY)
            if (length(summarY) == 1){ ## if there is only one level left
              oerr <- 0
              misY <- factor(rep(names(summarY), length(misi)))
            } else {
              RF <- local_linear_forest(
                X = obsX,
                Y = obsY,
                sample.fraction = 0.5, mtry = NULL,
                num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                alpha = NULL, imbalance.penalty = NULL,
                compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                samples_per_cluster = NULL, tune.parameters = F,
                num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
              ## record out-of-bag error
              oerr <-0 # RF$err.rate[[num.trees,1]]
              #             }
              ## predict missing values in column varInd
              misY <- predict(RF, misX)$predictions
            }
          }
          list(varInd=varInd, misY=misY, oerr=oerr)
        }
        ## update the master copy of the data
        for (res in results) {
          misi <- NAloc[,res$varInd]
          ximp[misi, res$varInd] <- res$misY
          OOBerror[res$varInd] <- res$oerr
        }
      }
    } else { # if parallelize != "variables"
      for (s in 1:p) {
        varInd <- sort.j[s]
        if (noNAvar[[varInd]] != 0) {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- ximp[obsi, varInd]
          obsX <- ximp[obsi, seq(1, p)[-varInd]]
          misX <- ximp[misi, seq(1, p)[-varInd]]
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            if (parallelize == 'forests') {
              xntree <- NULL
              RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                            .combine='combine', .multicombine=TRUE,
                            .packages='grf') %dopar% {
                              local_linear_forest(
                                X = obsX,
                                Y = obsY,
                                sample.fraction = 0.5, mtry = NULL,
                                num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                                honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                                alpha = NULL, imbalance.penalty = NULL,
                                compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                                samples_per_cluster = NULL, tune.parameters = F,
                                num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                            }
              ## record out-of-bag error
              OOBerror[varInd] <- mean((predict(RF) - RF$y) ^ 2, na.rm=TRUE)
              #               OOBerror[varInd] <- RF$mse[ntree]
            } else {
              RF <- local_linear_forest(
                X = obsX,
                Y = obsY,
                sample.fraction = 0.5, mtry = NULL,
                num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                alpha = NULL, imbalance.penalty = NULL,
                compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                samples_per_cluster = NULL, tune.parameters = F,
                num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
              ## record out-of-bag error
              OOBerror[varInd] <-0 # RF$mse[num.trees]
            }
            misY <- predict(RF, misX)$predictions
          } else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              misY <- factor(rep(names(summarY), sum(misi)))
            } else {
              if (parallelize == 'forests') {
                RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                              .combine='combine', .multicombine=TRUE,
                              .packages='grf') %dopar% {
                                local_linear_forest(
                                  X = obsX,
                                  Y = obsY,
                                  sample.fraction = 0.5, mtry = NULL,
                                  num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                                  honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                                  alpha = NULL, imbalance.penalty = NULL,
                                  compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                                  samples_per_cluster = NULL, tune.parameters = F,
                                  num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                              }
                ## record out-of-bag error
                ne <- as.integer(predict(RF)) != as.integer(RF$y)
                ne <- ne[! is.na(ne)]
                OOBerror[varInd] <- sum(ne) / length(ne)
              } else {
                RF <- local_linear_forest(
                  X = obsX,
                  Y = obsY,
                  sample.fraction = 0.5, mtry = NULL,
                  num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                  honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                  alpha = NULL, imbalance.penalty = NULL,
                  compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                  samples_per_cluster = NULL, tune.parameters = F,
                  num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                ## record out-of-bag error
                OOBerror[varInd] <-0 # RF$err.rate[[num.trees, 1]]
              }
              ## predict missing parts of Y
              misY <- predict(RF, misX)$predictions
            }
          }
          ximp[misi, varInd] <- misY
        }
      }
    }
    cat('done!\n')
    
    iter <- iter+1
    Ximp[[iter]] <- ximp
    
    t.co2 <- 1
    ## check the difference between iteration steps
    for (t.type in names(convNew)){
      t.ind <- which(varType == t.type)
      if (t.type == "numeric"){
        convNew[t.co2] <- sum((ximp[,t.ind]-ximp.old[,t.ind])^2)/sum(ximp[,t.ind]^2)
      } else {
        dist <- sum(as.character(as.matrix(ximp[,t.ind])) != as.character(as.matrix(ximp.old[,t.ind])))
        convNew[t.co2] <- dist / (n * sum(varType == 'factor'))
      }
      t.co2 <- t.co2 + 1
    }
    
    ## compute estimated imputation error
    if (!variablewise){
      NRMSE <- sqrt(mean(OOBerror[varType=='numeric'])/
                      var(as.vector(as.matrix(xmis[,varType=='numeric'])),
                          na.rm = TRUE))
      PFC <- mean(OOBerror[varType=='factor'])
      if (k==1){
        if (unique(varType)=='numeric'){
          OOBerr <- NRMSE
          names(OOBerr) <- 'NRMSE'
        } else {
          OOBerr <- PFC
          names(OOBerr) <- 'PFC'
        }
      } else {
        OOBerr <- c(NRMSE, PFC)
        names(OOBerr) <- c('NRMSE', 'PFC')
      }
    } else {
      OOBerr <- OOBerror
      names(OOBerr)[varType=='numeric'] <- 'MSE'
      names(OOBerr)[varType=='factor'] <- 'PFC'
    }
    
    if (any(!is.na(xtrue))){
      err <- suppressWarnings(mixError(ximp, xmis, xtrue))
    }
    
    ## return status output, if desired
    if (verbose){
      delta.start <- proc.time() - t.start
      if (any(!is.na(xtrue))){
        cat("    error(s):", err, "\n")
      }
      cat("    estimated error(s):", OOBerr, "\n")
      cat("    difference(s):", convNew, "\n")
      cat("    time:", delta.start[3], "seconds\n\n")
    }
  }#end while((convNew<convOld)&(iter<maxiter)){
  
  ## produce output w.r.t. stopping rule
  if (iter == maxiter){
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr)
    } else {
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, error = err)
    }
  } else {
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld)
    } else {
      out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld,
                  error = suppressWarnings(mixError(Ximp[[iter-1]], xmis, xtrue)))
    }
  }
  class(out) <- 'missForest'
  return(out)
}




missForest2 <- function(xmis, maxiter = 10, ntree = 100, variablewise = FALSE,
                        decreasing = FALSE, verbose = FALSE,
                        mtry = floor(sqrt(ncol(xmis))), replace = TRUE,
                        classwt = NULL, cutoff = NULL, strata = NULL,
                        sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                        xtrue = NA, parallelize = c('no', 'variables', 'forests'))
{ ## ----------------------------------------------------------------------
  ## Arguments:
  ## xmis         = data matrix with missing values
  ## maxiter      = stop after how many iterations (default = 10)
  ## ntree        = how many trees are grown in the forest (default = 100)
  ## variablewise = (boolean) return OOB errors for each variable separately
  ## decreasing   = (boolean) if TRUE the columns are sorted with decreasing
  ##                amount of missing values
  ## verbose      = (boolean) if TRUE then missForest returns error estimates,
  ##                runtime and if available true error during iterations
  ## mtry         = how many variables should be tried randomly at each node
  ## replace      = (boolean) if TRUE bootstrap sampling (with replacements)
  ##                is performed, else subsampling (without replacements)
  ## classwt      = list of priors of the classes in the categorical variables
  ## cutoff       = list of class cutoffs for each categorical variable
  ## strata       = list of (factor) variables used for stratified sampling
  ## sampsize     = list of size(s) of sample to draw
  ## nodesize     = minimum size of terminal nodes, vector of length 2, with
  ##                number for continuous variables in the first entry and
  ##                number for categorical variables in the second entry
  ## maxnodes     = maximum number of terminal nodes for individual trees
  ## xtrue        = complete data matrix
  ##
  ## ----------------------------------------------------------------------
  ## Author: Daniel Stekhoven, stekhoven@stat.math.ethz.ch
  
  ## stop in case of wrong inputs passed to randomForest
  n <- nrow(xmis)
  
  p <- ncol(xmis)
  if (!is.null(classwt))
    stopifnot(length(classwt) == p, typeof(classwt) == 'list')
  if (!is.null(cutoff))
    stopifnot(length(cutoff) == p, typeof(cutoff) == 'list')
  if (!is.null(strata))
    stopifnot(length(strata) == p, typeof(strata) == 'list')
  if (!is.null(nodesize))
    stopifnot(length(nodesize) == 2)
  
  ## remove completely missing variables
  if (any(apply(is.na(xmis), 2, sum) == n)){
    indCmis <- which(apply(is.na(xmis), 2, sum) == n)
    xmis <- xmis[,-indCmis]
    p <- ncol(xmis)
    cat('  removed variable(s)', indCmis,
        'due to the missingness of all entries\n')
  } 
  
  ## return feedback on parallelization setup
  parallelize <- match.arg(parallelize)
  if (parallelize %in% c('variables', 'forests')) {
    if (getDoParWorkers() == 1) {
      stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
    } else if (verbose) {
      if (parallelize == 'variables') {
        cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
      } else {
        cat("  parallelizing computation of the random forest model objects\n")
      }
    }
    if (getDoParWorkers() > p){
      stop('The number of parallel cores should not exceed the number of variables (p=', p, ")")
    }
  }
  
  ## perform initial S.W.A.G. on xmis (mean imputation)
  
  ximp <- xmis
  xAttrib <- lapply(xmis, attributes)
  
  p=m
  res= FPCA(mlisty, mlistx,
            list(dataType='Sparse', error=FALSE,nRegGrid = m, kernel='epan', verbose=TRUE))
  
  X.pace=fitted.values(res)
  varType <- character(p)
  for (t.co in 1:p){
    if (is.null(xAttrib[[t.co]])){
      varType[t.co] <- 'numeric'
      ximp[is.na(xmis[,t.co]),t.co] <- X.pace[is.na(xmis[,t.co]),t.co]
    } else {
      varType[t.co] <- 'factor'
      ## take the level which is more 'likely' (majority vote)
      max.level <- max(table(ximp[,t.co]))
      ## if there are several classes which are major, sample one at random
      class.assign <- sample(names(which(max.level == summary(ximp[,t.co]))), 1)
      ## it shouldn't be the NA class
      if (class.assign != "NA's"){
        ximp[is.na(xmis[,t.co]),t.co] <- class.assign
      } else {
        while (class.assign == "NA's"){
          class.assign <- sample(names(which(max.level ==
                                               summary(ximp[,t.co]))), 1)
        }
        ximp[is.na(xmis[,t.co]),t.co] <- class.assign
      }
    }
  }
  
  
  ## extract missingness pattern
  NAloc <- is.na(xmis)            # where are missings
  noNAvar <- apply(NAloc, 2, sum) # how many are missing in the vars
  sort.j <- order(noNAvar)        # indices of increasing amount of NA in vars
  if (decreasing)
    sort.j <- rev(sort.j)
  sort.noNAvar <- noNAvar[sort.j]
  
  ## compute a list of column indices for variable parallelization
  nzsort.j <- sort.j[sort.noNAvar > 0]
  if (parallelize == 'variables') {
    '%cols%' <- get('%dopar%')
    idxList <- as.list(isplitVector(nzsort.j, chunkSize=getDoParWorkers()))
  } 
  #   else {
  #     ## force column loop to be sequential
  #     '%cols%' <- get('%do%')
  #     idxList <- nzsort.j
  #   }
  
  ## output
  Ximp <- vector('list', maxiter)
  
  ## initialize parameters of interest
  iter <- 0
  k <- length(unique(varType))
  convNew <- rep(0, k)
  convOld <- rep(Inf, k)
  OOBerror <- numeric(p)
  names(OOBerror) <- varType
  
  ## setup convergence variables w.r.t. variable types
  if (k == 1){
    if (unique(varType) == 'numeric'){
      names(convNew) <- c('numeric')
    } else {
      names(convNew) <- c('factor')
    }
    convergence <- c()
    OOBerr <- numeric(1)
  } else {
    names(convNew) <- c('numeric', 'factor')
    convergence <- matrix(NA, ncol = 2)
    OOBerr <- numeric(2)
  }
  
  ## function to yield the stopping criterion in the following 'while' loop
  stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
    k <- length(unique(varType))
    if (k == 1){
      (convNew < convOld) & (iter < maxiter)
    } else {
      ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)
    }
  }
  
  ## iterate missForest
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)){
    if (iter != 0){
      convOld <- convNew
      OOBerrOld <- OOBerr
    }
    cat("  missForest iteration", iter+1, "in progress...")
    t.start <- proc.time()
    ximp.old <- ximp
    
    if (parallelize=="variables"){
      for (idx in idxList) {
        results <- foreach(varInd=idx, .packages='grf') %cols% {
          obsi <- !NAloc[,varInd] # which i's are observed
          misi <- NAloc[,varInd] # which i's are missing
          obsY <- ximp[obsi, varInd] # training response
          obsX <- ximp[obsi, seq(1, p)[-varInd]] # training variables
          misX <- ximp[misi, seq(1, p)[-varInd]] # prediction variables
          typeY <- varType[varInd]
          if (typeY == 'numeric'){
            RF <- local_linear_forest(
              X = obsX,
              Y = obsY,
              sample.fraction = 0.5, mtry = NULL,
              num.trees = 2000, num.threads = NULL, min.node.size = NULL,
              honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
              alpha = NULL, imbalance.penalty = NULL,
              compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
              samples_per_cluster = NULL, tune.parameters = FALSE,
              num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
            ## record out-of-bag error
            oerr <- 0 #RF$mse[num.trees]
            #           }
            ## predict missing values in column varInd
            misY <- predict(RF, misX)$predictions
          } else { # if Y is categorical          
            obsY <- factor(obsY) ## remove empty classes
            summarY <- summary(obsY)
            if (length(summarY) == 1){ ## if there is only one level left
              oerr <- 0
              misY <- factor(rep(names(summarY), length(misi)))
            } else {
              RF <- local_linear_forest(
                X = obsX,
                Y = obsY,
                sample.fraction = 0.5, mtry = NULL,
                num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                alpha = NULL, imbalance.penalty = NULL,
                compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                samples_per_cluster = NULL, tune.parameters = FALSE,
                num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
              ## record out-of-bag error
              oerr <-0 # RF$err.rate[[num.trees,1]]
              #             }
              ## predict missing values in column varInd
              misY <- predict(RF, misX)$predictions
            }
          }
          list(varInd=varInd, misY=misY, oerr=oerr)
        }
        ## update the master copy of the data
        for (res in results) {
          misi <- NAloc[,res$varInd]
          ximp[misi, res$varInd] <- res$misY
          OOBerror[res$varInd] <- res$oerr
        }
      }
    } else { # if parallelize != "variables"
      for (s in 1:p) {
        varInd <- sort.j[s]
        if (noNAvar[[varInd]] != 0) {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- ximp[obsi, varInd]
          obsX <- ximp[obsi, seq(1, p)[-varInd]]
          misX <- ximp[misi, seq(1, p)[-varInd]]
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            if (parallelize == 'forests') {
              xntree <- NULL
              RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                            .combine='combine', .multicombine=TRUE,
                            .packages='grf') %dopar% {
                              local_linear_forest(
                                X = obsX,
                                Y = obsY,
                                sample.fraction = 0.5, mtry = NULL,
                                num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                                honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                                alpha = NULL, imbalance.penalty = NULL,
                                compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                                samples_per_cluster = NULL, tune.parameters = FALSE,
                                num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                            }
              ## record out-of-bag error
              OOBerror[varInd] <- mean((predict(RF) - RF$y) ^ 2, na.rm=TRUE)
              #               OOBerror[varInd] <- RF$mse[ntree]
            } else {
              RF <- local_linear_forest(
                X = obsX,
                Y = obsY,
                sample.fraction = 0.5, mtry = NULL,
                num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                alpha = NULL, imbalance.penalty = NULL,
                compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                samples_per_cluster = NULL, tune.parameters = FALSE,
                num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
              ## record out-of-bag error
              OOBerror[varInd] <-0 # RF$mse[num.trees]
            }
            misY <- predict(RF, misX)$predictions
          } else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              misY <- factor(rep(names(summarY), sum(misi)))
            } else {
              if (parallelize == 'forests') {
                RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                              .combine='combine', .multicombine=TRUE,
                              .packages='grf') %dopar% {
                                local_linear_forest(
                                  X = obsX,
                                  Y = obsY,
                                  sample.fraction = 0.5, mtry = NULL,
                                  num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                                  honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                                  alpha = NULL, imbalance.penalty = NULL,
                                  compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                                  samples_per_cluster = NULL, tune.parameters = FALSE,
                                  num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                              }
                ## record out-of-bag error
                ne <- as.integer(predict(RF)) != as.integer(RF$y)
                ne <- ne[! is.na(ne)]
                OOBerror[varInd] <- sum(ne) / length(ne)
              } else {
                RF <- local_linear_forest(
                  X = obsX,
                  Y = obsY,
                  sample.fraction = 0.5, mtry = NULL,
                  num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                  honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                  alpha = NULL, imbalance.penalty = NULL,
                  compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                  samples_per_cluster = NULL, tune.parameters = FALSE,
                  num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                ## record out-of-bag error
                OOBerror[varInd] <-0 # RF$err.rate[[num.trees, 1]]
              }
              ## predict missing parts of Y
              misY <- predict(RF, misX)$predictions
            }
          }
          ximp[misi, varInd] <- misY
        }
      }
    }
    cat('done!\n')
    
    iter <- iter+1
    Ximp[[iter]] <- ximp
    
    t.co2 <- 1
    ## check the difference between iteration steps
    for (t.type in names(convNew)){
      t.ind <- which(varType == t.type)
      if (t.type == "numeric"){
        convNew[t.co2] <- sum((ximp[,t.ind]-ximp.old[,t.ind])^2)/sum(ximp[,t.ind]^2)
      } else {
        dist <- sum(as.character(as.matrix(ximp[,t.ind])) != as.character(as.matrix(ximp.old[,t.ind])))
        convNew[t.co2] <- dist / (n * sum(varType == 'factor'))
      }
      t.co2 <- t.co2 + 1
    }
    
    ## compute estimated imputation error
    if (!variablewise){
      NRMSE <- sqrt(mean(OOBerror[varType=='numeric'])/
                      var(as.vector(as.matrix(xmis[,varType=='numeric'])),
                          na.rm = TRUE))
      PFC <- mean(OOBerror[varType=='factor'])
      if (k==1){
        if (unique(varType)=='numeric'){
          OOBerr <- NRMSE
          names(OOBerr) <- 'NRMSE'
        } else {
          OOBerr <- PFC
          names(OOBerr) <- 'PFC'
        }
      } else {
        OOBerr <- c(NRMSE, PFC)
        names(OOBerr) <- c('NRMSE', 'PFC')
      }
    } else {
      OOBerr <- OOBerror
      names(OOBerr)[varType=='numeric'] <- 'MSE'
      names(OOBerr)[varType=='factor'] <- 'PFC'
    }
    
    if (any(!is.na(xtrue))){
      err <- suppressWarnings(mixError(ximp, xmis, xtrue))
    }
    
    ## return status output, if desired
    if (verbose){
      delta.start <- proc.time() - t.start
      if (any(!is.na(xtrue))){
        cat("    error(s):", err, "\n")
      }
      cat("    estimated error(s):", OOBerr, "\n")
      cat("    difference(s):", convNew, "\n")
      cat("    time:", delta.start[3], "seconds\n\n")
    }
  }#end while((convNew<convOld)&(iter<maxiter)){
  
  ## produce output w.r.t. stopping rule
  if (iter == maxiter){
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr)
    } else {
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, error = err)
    }
  } else {
    if (any(is.na(xtrue))){
      out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld)
    } else {
      out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld,
                  error = suppressWarnings(mixError(Ximp[[iter-1]], xmis, xtrue)))
    }
  }
  class(out) <- 'missForest'
  return(out)
}








allmodelimputeb=function(n,m,o,seeds,buck){
  
  
  pacemissForestb <- function(xmis, maxiter = 10, ntree = 100, variablewise = FALSE,
                              decreasing = FALSE, verbose = FALSE,
                              mtry = floor(sqrt(ncol(xmis))), replace = TRUE,
                              classwt = NULL, cutoff = NULL, strata = NULL,
                              sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                              xtrue = NA, parallelize = c('no', 'variables', 'forests'))
  { 
    n <- nrow(xmis)
    p <- ncol(xmis)
    if (!is.null(classwt))
      stopifnot(length(classwt) == p, typeof(classwt) == 'list')
    if (!is.null(cutoff))
      stopifnot(length(cutoff) == p, typeof(cutoff) == 'list')
    if (!is.null(strata))
      stopifnot(length(strata) == p, typeof(strata) == 'list')
    if (!is.null(nodesize))
      stopifnot(length(nodesize) == 2)
    
    ## remove completely missing variables
    if (any(apply(is.na(xmis), 2, sum) == n)){
      indCmis <- which(apply(is.na(xmis), 2, sum) == n)
      xmis <- xmis[,-indCmis]
      p <- ncol(xmis)
      cat('  removed variable(s)', indCmis,
          'due to the missingness of all entries\n')
    } 
    
    ## return feedback on parallelization setup
    parallelize <- match.arg(parallelize)
    if (parallelize %in% c('variables', 'forests')) {
      if (getDoParWorkers() == 1) {
        stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
      } else if (verbose) {
        if (parallelize == 'variables') {
          cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
        } else {
          cat("  parallelizing computation of the random forest model objects\n")
        }
      }
      if (getDoParWorkers() > p){
        stop('The number of parallel cores should not exceed the number of variables (p=', p, ")")
      }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    ## perform initial S.W.A.G. on xmis (mean imputation)
    
    ximp <- xmis
    xAttrib <- lapply(xmis, attributes)
    
    p=buck
    res= FPCA(mlisty, mlistx,
              list(dataType='Sparse', error=FALSE,nRegGrid = buck, kernel='epan', verbose=TRUE,userBwCov = 2))
    
    X.pace=fitted.values(res)
    varType <- character(p)
    for (t.co in 1:p){
      if (is.null(xAttrib[[t.co]])){
        varType[t.co] <- 'numeric'
        ximp[is.na(xmis[,t.co]),t.co] <- X.pace[is.na(xmis[,t.co]),t.co]
      } else {
        varType[t.co] <- 'factor'
        ## take the level which is more 'likely' (majority vote)
        max.level <- max(table(ximp[,t.co]))
        ## if there are several classes which are major, sample one at random
        class.assign <- sample(names(which(max.level == summary(ximp[,t.co]))), 1)
        ## it shouldn't be the NA class
        if (class.assign != "NA's"){
          ximp[is.na(xmis[,t.co]),t.co] <- class.assign
        } else {
          while (class.assign == "NA's"){
            class.assign <- sample(names(which(max.level ==
                                                 summary(ximp[,t.co]))), 1)
          }
          ximp[is.na(xmis[,t.co]),t.co] <- class.assign
        }
      }
    }
    
    ## extract missingness pattern
    NAloc <- is.na(xmis)            # where are missings
    noNAvar <- apply(NAloc, 2, sum) # how many are missing in the vars
    sort.j <- order(noNAvar)        # indices of increasing amount of NA in vars
    if (decreasing)
      sort.j <- rev(sort.j)
    sort.noNAvar <- noNAvar[sort.j]
    
    ## compute a list of column indices for variable parallelization
    nzsort.j <- sort.j[sort.noNAvar > 0]
    if (parallelize == 'variables') {
      '%cols%' <- get('%dopar%')
      idxList <- as.list(isplitVector(nzsort.j, chunkSize=getDoParWorkers()))
    } 
    #   else {
    #     ## force column loop to be sequential
    #     '%cols%' <- get('%do%')
    #     idxList <- nzsort.j
    #   }
    
    ## output
    Ximp <- vector('list', maxiter)
    
    ## initialize parameters of interest
    iter <- 0
    k <- length(unique(varType))
    convNew <- rep(0, k)
    convOld <- rep(Inf, k)
    OOBerror <- numeric(p)
    names(OOBerror) <- varType
    
    ## setup convergence variables w.r.t. variable types
    if (k == 1){
      if (unique(varType) == 'numeric'){
        names(convNew) <- c('numeric')
      } else {
        names(convNew) <- c('factor')
      }
      convergence <- c()
      OOBerr <- numeric(1)
    } else {
      names(convNew) <- c('numeric', 'factor')
      convergence <- matrix(NA, ncol = 2)
      OOBerr <- numeric(2)
    }
    
    ## function to yield the stopping criterion in the following 'while' loop
    stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
      k <- length(unique(varType))
      if (k == 1){
        (convNew < convOld) & (iter < maxiter)
      } else {
        ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)
      }
    }
    
    ## iterate missForest
    while (stopCriterion(varType, convNew, convOld, iter, maxiter)){
      if (iter != 0){
        convOld <- convNew
        OOBerrOld <- OOBerr
      }
      cat("  missForest iteration", iter+1, "in progress...")
      t.start <- proc.time()
      ximp.old <- ximp
      
      if (parallelize=="variables"){
        for (idx in idxList) {
          results <- foreach(varInd=idx, .packages='randomForest') %cols% {
            obsi <- !NAloc[,varInd] # which i's are observed
            misi <- NAloc[,varInd] # which i's are missing
            obsY <- ximp[obsi, varInd] # training response
            obsX <- ximp[obsi, seq(1, p)[-varInd]] # training variables
            misX <- ximp[misi, seq(1, p)[-varInd]] # prediction variables
            typeY <- varType[varInd]
            if (typeY == 'numeric'){
              RF <- randomForest(
                x = obsX,
                y = obsY,
                ntree = ntree,
                mtry = mtry,
                replace = replace,
                sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                  if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
              ## record out-of-bag error
              oerr <- RF$mse[ntree]
              #           }
              ## predict missing values in column varInd
              misY <- predict(RF, misX)
            } else { # if Y is categorical          
              obsY <- factor(obsY) ## remove empty classes
              summarY <- summary(obsY)
              if (length(summarY) == 1){ ## if there is only one level left
                oerr <- 0
                misY <- factor(rep(names(summarY), length(misi)))
              } else {
                RF <- randomForest(
                  x = obsX,
                  y = obsY,
                  ntree = ntree,
                  mtry = mtry,
                  replace = replace,
                  classwt = if (!is.null(classwt)) classwt[[varInd]] else
                    rep(1, nlevels(obsY)),
                  cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else
                    rep(1/nlevels(obsY), nlevels(obsY)),
                  strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                  sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                    if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                  nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                  maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                ## record out-of-bag error
                oerr <- RF$err.rate[[ntree,1]]
                #             }
                ## predict missing values in column varInd
                misY <- predict(RF, misX)
              }
            }
            list(varInd=varInd, misY=misY, oerr=oerr)
          }
          ## update the master copy of the data
          for (res in results) {
            misi <- NAloc[,res$varInd]
            ximp[misi, res$varInd] <- res$misY
            OOBerror[res$varInd] <- res$oerr
          }
        }
      } else { # if parallelize != "variables"
        for (s in 1:p) {
          varInd <- sort.j[s]
          if (noNAvar[[varInd]] != 0) {
            obsi <- !NAloc[, varInd]
            misi <- NAloc[, varInd]
            obsY <- ximp[obsi, varInd]
            obsX <- ximp[obsi, seq(1, p)[-varInd]]
            misX <- ximp[misi, seq(1, p)[-varInd]]
            typeY <- varType[varInd]
            if (typeY == "numeric") {
              if (parallelize == 'forests') {
                xntree <- NULL
                RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                              .combine='combine', .multicombine=TRUE,
                              .packages='randomForest') %dopar% {
                                randomForest( x = obsX,
                                              y = obsY,
                                              ntree = xntree,
                                              mtry = mtry,
                                              replace = replace,
                                              sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                                if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                                              nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                                              maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                              }
                ## record out-of-bag error
                OOBerror[varInd] <- mean((predict(RF) - RF$y) ^ 2, na.rm=TRUE)
                #               OOBerror[varInd] <- RF$mse[ntree]
              } else {
                RF <- randomForest( x = obsX,
                                    y = obsY,
                                    ntree = ntree,
                                    mtry = mtry,
                                    replace = replace,
                                    sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                      if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                                    nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                                    maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                ## record out-of-bag error
                OOBerror[varInd] <- RF$mse[ntree]
              }
              misY <- predict(RF, misX)
            } else {
              obsY <- factor(obsY)
              summarY <- summary(obsY)
              if (length(summarY) == 1) {
                misY <- factor(rep(names(summarY), sum(misi)))
              } else {
                if (parallelize == 'forests') {
                  RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                                .combine='combine', .multicombine=TRUE,
                                .packages='randomForest') %dopar% {
                                  randomForest(
                                    x = obsX,
                                    y = obsY,
                                    ntree = xntree,
                                    mtry = mtry,
                                    replace = replace,
                                    classwt = if (!is.null(classwt)) classwt[[varInd]] else
                                      rep(1, nlevels(obsY)),
                                    cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else
                                      rep(1/nlevels(obsY), nlevels(obsY)),
                                    strata = if (!is.null(strata)) strata[[varInd]] else obsY,
                                    sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else
                                      if (replace) nrow(obsX) else ceiling(0.632*nrow(obsX)),
                                    nodesize = if (!is.null(nodesize)) nodesize[2] else 5,
                                    maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                                }
                  ## record out-of-bag error
                  ne <- as.integer(predict(RF)) != as.integer(RF$y)
                  ne <- ne[! is.na(ne)]
                  OOBerror[varInd] <- sum(ne) / length(ne)
                } else {
                  RF <- randomForest(x = obsX, 
                                     y = obsY, 
                                     ntree = ntree, 
                                     mtry = mtry, 
                                     replace = replace, 
                                     classwt = if (!is.null(classwt)) classwt[[varInd]] else 
                                       rep(1, nlevels(obsY)),
                                     cutoff = if (!is.null(cutoff)) cutoff[[varInd]] else 
                                       rep(1/nlevels(obsY), nlevels(obsY)),
                                     strata = if (!is.null(strata)) strata[[varInd]] else obsY, 
                                     sampsize = if (!is.null(sampsize)) sampsize[[varInd]] else 
                                       if (replace) nrow(obsX) else ceiling(0.632 * nrow(obsX)), 
                                     nodesize = if (!is.null(nodesize)) nodesize[2] else 5, 
                                     maxnodes = if (!is.null(maxnodes)) maxnodes else NULL)
                  ## record out-of-bag error
                  OOBerror[varInd] <- RF$err.rate[[ntree, 1]]
                }
                ## predict missing parts of Y
                misY <- predict(RF, misX)
              }
            }
            ximp[misi, varInd] <- misY
          }
        }
      }
      cat('done!\n')
      
      iter <- iter+1
      Ximp[[iter]] <- ximp
      
      t.co2 <- 1
      ## check the difference between iteration steps
      for (t.type in names(convNew)){
        t.ind <- which(varType == t.type)
        if (t.type == "numeric"){
          convNew[t.co2] <- sum((ximp[,t.ind]-ximp.old[,t.ind])^2)/sum(ximp[,t.ind]^2)
        } else {
          dist <- sum(as.character(as.matrix(ximp[,t.ind])) != as.character(as.matrix(ximp.old[,t.ind])))
          convNew[t.co2] <- dist / (n * sum(varType == 'factor'))
        }
        t.co2 <- t.co2 + 1
      }
      
      ## compute estimated imputation error
      if (!variablewise){
        NRMSE <- sqrt(mean(OOBerror[varType=='numeric'])/
                        var(as.vector(as.matrix(xmis[,varType=='numeric'])),
                            na.rm = TRUE))
        PFC <- mean(OOBerror[varType=='factor'])
        if (k==1){
          if (unique(varType)=='numeric'){
            OOBerr <- NRMSE
            names(OOBerr) <- 'NRMSE'
          } else {
            OOBerr <- PFC
            names(OOBerr) <- 'PFC'
          }
        } else {
          OOBerr <- c(NRMSE, PFC)
          names(OOBerr) <- c('NRMSE', 'PFC')
        }
      } else {
        OOBerr <- OOBerror
        names(OOBerr)[varType=='numeric'] <- 'MSE'
        names(OOBerr)[varType=='factor'] <- 'PFC'
      }
      
      if (any(!is.na(xtrue))){
        err <- suppressWarnings(mixError(ximp, xmis, xtrue))
      }
      
      ## return status output, if desired
      if (verbose){
        delta.start <- proc.time() - t.start
        if (any(!is.na(xtrue))){
          cat("    error(s):", err, "\n")
        }
        cat("    estimated error(s):", OOBerr, "\n")
        cat("    difference(s):", convNew, "\n")
        cat("    time:", delta.start[3], "seconds\n\n")
      }
    }#end while((convNew<convOld)&(iter<maxiter)){
    
    ## produce output w.r.t. stopping rule
    if (iter == maxiter){
      if (any(is.na(xtrue))){
        out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr)
      } else {
        out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, error = err)
      }
    } else {
      if (any(is.na(xtrue))){
        out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld)
      } else {
        out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld,
                    error = suppressWarnings(mixError(Ximp[[iter-1]], xmis, xtrue)))
      }
    }
    class(out) <- 'missForest'
    return(out)
  }
  
  
  missForest2b <- function(xmis, maxiter = 10, ntree = 100, variablewise = FALSE,
                           decreasing = FALSE, verbose = FALSE,
                           mtry = floor(sqrt(ncol(xmis))), replace = TRUE,
                           classwt = NULL, cutoff = NULL, strata = NULL,
                           sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                           xtrue = NA, parallelize = c('no', 'variables', 'forests'))
  { ## ----------------------------------------------------------------------
    ## Arguments:
    ## xmis         = data matrix with missing values
    ## maxiter      = stop after how many iterations (default = 10)
    ## ntree        = how many trees are grown in the forest (default = 100)
    ## variablewise = (boolean) return OOB errors for each variable separately
    ## decreasing   = (boolean) if TRUE the columns are sorted with decreasing
    ##                amount of missing values
    ## verbose      = (boolean) if TRUE then missForest returns error estimates,
    ##                runtime and if available true error during iterations
    ## mtry         = how many variables should be tried randomly at each node
    ## replace      = (boolean) if TRUE bootstrap sampling (with replacements)
    ##                is performed, else subsampling (without replacements)
    ## classwt      = list of priors of the classes in the categorical variables
    ## cutoff       = list of class cutoffs for each categorical variable
    ## strata       = list of (factor) variables used for stratified sampling
    ## sampsize     = list of size(s) of sample to draw
    ## nodesize     = minimum size of terminal nodes, vector of length 2, with
    ##                number for continuous variables in the first entry and
    ##                number for categorical variables in the second entry
    ## maxnodes     = maximum number of terminal nodes for individual trees
    ## xtrue        = complete data matrix
    ##
    ## ----------------------------------------------------------------------
    ## Author: Daniel Stekhoven, stekhoven@stat.math.ethz.ch
    
    ## stop in case of wrong inputs passed to randomForest
    n <- nrow(xmis)
    
    p <- ncol(xmis)
    if (!is.null(classwt))
      stopifnot(length(classwt) == p, typeof(classwt) == 'list')
    if (!is.null(cutoff))
      stopifnot(length(cutoff) == p, typeof(cutoff) == 'list')
    if (!is.null(strata))
      stopifnot(length(strata) == p, typeof(strata) == 'list')
    if (!is.null(nodesize))
      stopifnot(length(nodesize) == 2)
    
    ## remove completely missing variables
    if (any(apply(is.na(xmis), 2, sum) == n)){
      indCmis <- which(apply(is.na(xmis), 2, sum) == n)
      xmis <- xmis[,-indCmis]
      p <- ncol(xmis)
      cat('  removed variable(s)', indCmis,
          'due to the missingness of all entries\n')
    } 
    
    ## return feedback on parallelization setup
    parallelize <- match.arg(parallelize)
    if (parallelize %in% c('variables', 'forests')) {
      if (getDoParWorkers() == 1) {
        stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
      } else if (verbose) {
        if (parallelize == 'variables') {
          cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
        } else {
          cat("  parallelizing computation of the random forest model objects\n")
        }
      }
      if (getDoParWorkers() > p){
        stop('The number of parallel cores should not exceed the number of variables (p=', p, ")")
      }
    }
    
    ## perform initial S.W.A.G. on xmis (mean imputation)
    
    ximp <- xmis
    xAttrib <- lapply(xmis, attributes)
    
    p=buck
    res= FPCA(mlisty, mlistx,
              list(dataType='Sparse', error=FALSE,nRegGrid = buck, kernel='epan', verbose=TRUE,userBwCov = 2))
    
    X.pace=fitted.values(res)
    varType <- character(p)
    for (t.co in 1:p){
      if (is.null(xAttrib[[t.co]])){
        varType[t.co] <- 'numeric'
        ximp[is.na(xmis[,t.co]),t.co] <- X.pace[is.na(xmis[,t.co]),t.co]
      } else {
        varType[t.co] <- 'factor'
        ## take the level which is more 'likely' (majority vote)
        max.level <- max(table(ximp[,t.co]))
        ## if there are several classes which are major, sample one at random
        class.assign <- sample(names(which(max.level == summary(ximp[,t.co]))), 1)
        ## it shouldn't be the NA class
        if (class.assign != "NA's"){
          ximp[is.na(xmis[,t.co]),t.co] <- class.assign
        } else {
          while (class.assign == "NA's"){
            class.assign <- sample(names(which(max.level ==
                                                 summary(ximp[,t.co]))), 1)
          }
          ximp[is.na(xmis[,t.co]),t.co] <- class.assign
        }
      }
    }
    
    
    ## extract missingness pattern
    NAloc <- is.na(xmis)            # where are missings
    noNAvar <- apply(NAloc, 2, sum) # how many are missing in the vars
    sort.j <- order(noNAvar)        # indices of increasing amount of NA in vars
    if (decreasing)
      sort.j <- rev(sort.j)
    sort.noNAvar <- noNAvar[sort.j]
    
    ## compute a list of column indices for variable parallelization
    nzsort.j <- sort.j[sort.noNAvar > 0]
    if (parallelize == 'variables') {
      '%cols%' <- get('%dopar%')
      idxList <- as.list(isplitVector(nzsort.j, chunkSize=getDoParWorkers()))
    } 
    #   else {
    #     ## force column loop to be sequential
    #     '%cols%' <- get('%do%')
    #     idxList <- nzsort.j
    #   }
    
    ## output
    Ximp <- vector('list', maxiter)
    
    ## initialize parameters of interest
    iter <- 0
    k <- length(unique(varType))
    convNew <- rep(0, k)
    convOld <- rep(Inf, k)
    OOBerror <- numeric(p)
    names(OOBerror) <- varType
    
    ## setup convergence variables w.r.t. variable types
    if (k == 1){
      if (unique(varType) == 'numeric'){
        names(convNew) <- c('numeric')
      } else {
        names(convNew) <- c('factor')
      }
      convergence <- c()
      OOBerr <- numeric(1)
    } else {
      names(convNew) <- c('numeric', 'factor')
      convergence <- matrix(NA, ncol = 2)
      OOBerr <- numeric(2)
    }
    
    ## function to yield the stopping criterion in the following 'while' loop
    stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
      k <- length(unique(varType))
      if (k == 1){
        (convNew < convOld) & (iter < maxiter)
      } else {
        ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)
      }
    }
    
    ## iterate missForest
    while (stopCriterion(varType, convNew, convOld, iter, maxiter)){
      if (iter != 0){
        convOld <- convNew
        OOBerrOld <- OOBerr
      }
      cat("  missForest iteration", iter+1, "in progress...")
      t.start <- proc.time()
      ximp.old <- ximp
      
      if (parallelize=="variables"){
        for (idx in idxList) {
          results <- foreach(varInd=idx, .packages='grf') %cols% {
            obsi <- !NAloc[,varInd] # which i's are observed
            misi <- NAloc[,varInd] # which i's are missing
            obsY <- ximp[obsi, varInd] # training response
            obsX <- ximp[obsi, seq(1, p)[-varInd]] # training variables
            misX <- ximp[misi, seq(1, p)[-varInd]] # prediction variables
            typeY <- varType[varInd]
            if (typeY == 'numeric'){
              RF <- local_linear_forest(
                X = obsX,
                Y = obsY,
                sample.fraction = 0.5, mtry = NULL,
                num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                alpha = NULL, imbalance.penalty = NULL,
                compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                samples_per_cluster = NULL, tune.parameters = FALSE,
                num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
              ## record out-of-bag error
              oerr <- 0 #RF$mse[num.trees]
              #           }
              ## predict missing values in column varInd
              misY <- predict(RF, misX)$predictions
            } else { # if Y is categorical          
              obsY <- factor(obsY) ## remove empty classes
              summarY <- summary(obsY)
              if (length(summarY) == 1){ ## if there is only one level left
                oerr <- 0
                misY <- factor(rep(names(summarY), length(misi)))
              } else {
                RF <- local_linear_forest(
                  X = obsX,
                  Y = obsY,
                  sample.fraction = 0.5, mtry = NULL,
                  num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                  honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                  alpha = NULL, imbalance.penalty = NULL,
                  compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                  samples_per_cluster = NULL, tune.parameters = FALSE,
                  num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                ## record out-of-bag error
                oerr <-0 # RF$err.rate[[num.trees,1]]
                #             }
                ## predict missing values in column varInd
                misY <- predict(RF, misX)$predictions
              }
            }
            list(varInd=varInd, misY=misY, oerr=oerr)
          }
          ## update the master copy of the data
          for (res in results) {
            misi <- NAloc[,res$varInd]
            ximp[misi, res$varInd] <- res$misY
            OOBerror[res$varInd] <- res$oerr
          }
        }
      } else { # if parallelize != "variables"
        for (s in 1:p) {
          varInd <- sort.j[s]
          if (noNAvar[[varInd]] != 0) {
            obsi <- !NAloc[, varInd]
            misi <- NAloc[, varInd]
            obsY <- ximp[obsi, varInd]
            obsX <- ximp[obsi, seq(1, p)[-varInd]]
            misX <- ximp[misi, seq(1, p)[-varInd]]
            typeY <- varType[varInd]
            if (typeY == "numeric") {
              if (parallelize == 'forests') {
                xntree <- NULL
                RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                              .combine='combine', .multicombine=TRUE,
                              .packages='grf') %dopar% {
                                local_linear_forest(
                                  X = obsX,
                                  Y = obsY,
                                  sample.fraction = 0.5, mtry = NULL,
                                  num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                                  honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                                  alpha = NULL, imbalance.penalty = NULL,
                                  compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                                  samples_per_cluster = NULL, tune.parameters = FALSE,
                                  num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                              }
                ## record out-of-bag error
                OOBerror[varInd] <- mean((predict(RF) - RF$y) ^ 2, na.rm=TRUE)
                #               OOBerror[varInd] <- RF$mse[ntree]
              } else {
                RF <- local_linear_forest(
                  X = obsX,
                  Y = obsY,
                  sample.fraction = 0.5, mtry = NULL,
                  num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                  honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                  alpha = NULL, imbalance.penalty = NULL,
                  compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                  samples_per_cluster = NULL, tune.parameters = FALSE,
                  num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                ## record out-of-bag error
                OOBerror[varInd] <-0 # RF$mse[num.trees]
              }
              misY <- predict(RF, misX)$predictions
            } else {
              obsY <- factor(obsY)
              summarY <- summary(obsY)
              if (length(summarY) == 1) {
                misY <- factor(rep(names(summarY), sum(misi)))
              } else {
                if (parallelize == 'forests') {
                  RF <- foreach(xntree=idiv(ntree, chunks=getDoParWorkers()),
                                .combine='combine', .multicombine=TRUE,
                                .packages='grf') %dopar% {
                                  local_linear_forest(
                                    X = obsX,
                                    Y = obsY,
                                    sample.fraction = 0.5, mtry = NULL,
                                    num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                                    honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                                    alpha = NULL, imbalance.penalty = NULL,
                                    compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                                    samples_per_cluster = NULL, tune.parameters = FALSE,
                                    num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                                }
                  ## record out-of-bag error
                  ne <- as.integer(predict(RF)) != as.integer(RF$y)
                  ne <- ne[! is.na(ne)]
                  OOBerror[varInd] <- sum(ne) / length(ne)
                } else {
                  RF <- local_linear_forest(
                    X = obsX,
                    Y = obsY,
                    sample.fraction = 0.5, mtry = NULL,
                    num.trees = 2000, num.threads = NULL, min.node.size = NULL,
                    honesty = TRUE, honesty.fraction = NULL, ci.group.size = 1,
                    alpha = NULL, imbalance.penalty = NULL,
                    compute.oob.predictions = FALSE, seed = NULL, clusters = NULL,
                    samples_per_cluster = NULL, tune.parameters = FALSE,
                    num.fit.trees = 10, num.fit.reps = 100, num.optimize.reps = 1000)
                  ## record out-of-bag error
                  OOBerror[varInd] <-0 # RF$err.rate[[num.trees, 1]]
                }
                ## predict missing parts of Y
                misY <- predict(RF, misX)$predictions
              }
            }
            ximp[misi, varInd] <- misY
          }
        }
      }
      cat('done!\n')
      
      iter <- iter+1
      Ximp[[iter]] <- ximp
      
      t.co2 <- 1
      ## check the difference between iteration steps
      for (t.type in names(convNew)){
        t.ind <- which(varType == t.type)
        if (t.type == "numeric"){
          convNew[t.co2] <- sum((ximp[,t.ind]-ximp.old[,t.ind])^2)/sum(ximp[,t.ind]^2)
        } else {
          dist <- sum(as.character(as.matrix(ximp[,t.ind])) != as.character(as.matrix(ximp.old[,t.ind])))
          convNew[t.co2] <- dist / (n * sum(varType == 'factor'))
        }
        t.co2 <- t.co2 + 1
      }
      
      ## compute estimated imputation error
      if (!variablewise){
        NRMSE <- sqrt(mean(OOBerror[varType=='numeric'])/
                        var(as.vector(as.matrix(xmis[,varType=='numeric'])),
                            na.rm = TRUE))
        PFC <- mean(OOBerror[varType=='factor'])
        if (k==1){
          if (unique(varType)=='numeric'){
            OOBerr <- NRMSE
            names(OOBerr) <- 'NRMSE'
          } else {
            OOBerr <- PFC
            names(OOBerr) <- 'PFC'
          }
        } else {
          OOBerr <- c(NRMSE, PFC)
          names(OOBerr) <- c('NRMSE', 'PFC')
        }
      } else {
        OOBerr <- OOBerror
        names(OOBerr)[varType=='numeric'] <- 'MSE'
        names(OOBerr)[varType=='factor'] <- 'PFC'
      }
      
      if (any(!is.na(xtrue))){
        err <- suppressWarnings(mixError(ximp, xmis, xtrue))
      }
      
      ## return status output, if desired
      if (verbose){
        delta.start <- proc.time() - t.start
        if (any(!is.na(xtrue))){
          cat("    error(s):", err, "\n")
        }
        cat("    estimated error(s):", OOBerr, "\n")
        cat("    difference(s):", convNew, "\n")
        cat("    time:", delta.start[3], "seconds\n\n")
      }
    }#end while((convNew<convOld)&(iter<maxiter)){
    
    ## produce output w.r.t. stopping rule
    if (iter == maxiter){
      if (any(is.na(xtrue))){
        out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr)
      } else {
        out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, error = err)
      }
    } else {
      if (any(is.na(xtrue))){
        out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld)
      } else {
        out <- list(ximp = Ximp[[iter-1]], OOBerror = OOBerrOld,
                    error = suppressWarnings(mixError(Ximp[[iter-1]], xmis, xtrue)))
      }
    }
    class(out) <- 'missForest'
    return(out)
  }
  
  
  
  x = seq(0,1,len=m)        #gird values
  
  cov<-function(t,s,sig2=1,rho=0.5){ # matern ftn nu = 5/2
    d <- abs(outer(t,s,"-"))
    tmp2 <- sig2*(1+sqrt(5)*d/rho + 5*d^2/(3*rho^2))*exp(-sqrt(5)*d/rho)}
  COV = cov(x, x)
  
  set.seed(seeds)
  data= mvrnorm(n, rep(0, length=length(x)), COV) #+ .3*mvrnorm(n, rep(0, length=length(x)), .1*diag(length(x))) #error
  
  data1=data.frame(data)
  
  dat= data.frame(x=x, t(data))
  dat= melt(dat, id="x")
  
  
  #plot for all the points of the data
  fig1= ggplot(dat,aes(x=x,y=value)) +
    geom_line(aes(group=variable)) +   theme_bw() +
    geom_line(aes(color=variable),show.legend = FALSE) + #REPLICATES
    
    scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
    xlab("input, x")+
    ggtitle("Plot of the data")
  fig1
  
  
  #---------------------------------------------------------------------------------------------
  #Imputation
  
  
  
  #generationg missingness
  set.seed(seeds)
  oo=m-floor(m*o)
  maa=c(rep(NA,m-oo),rep(1,oo))
  mb=matrix(NA,n,m)
  for(i in 1:n){
    mb[i,]=sample(maa)
  }
  mdata=mb*data  #MCAR
  
  mdata1=data.frame(mdata)
  
  mdat= data.frame(x=x, t(mdata1))
  mdat= melt(mdat, id="x")
  
  
  bucket=function(data,buck,timepoints,m){
    
    splitdf <- function(df, n) {
      indx <- matrix(seq_len(ncol(df)), ncol = n)
      lapply(seq_len(n), function(x) df[, indx[, x]])
    }
    
    a=splitdf((data),buck)
    
    aa=matrix(NA,n,buck)
    for(i in 1:buck) {
      # i-th element of `u1` squared into `i`-th position of `usq`
      aa[,i] <- rowMeans(a[[i]], na.rm=TRUE)
      
    }
    x1= seq_along(timepoints)
    d=split(timepoints,ceiling(x1/(m/buck)))
    dd1=d[[1]]
    dd2=d[[2]]
    
    
    xnew=cbind(mean(dd1),mean(dd2))
    
    xnew=matrix(NA,1,buck)
    for(i in 1:buck) {
      # i-th element of `u1` squared into `i`-th position of `usq`
      xnew[,i] <- mean(d[[i]])
      
    }
    return(list(aa,xnew))
  }
  
  
  bb=buck-2
  
  databucket=bucket(data=mdata[,3:m-1],buck=bb,timepoints=x[3:m-1],m=m-2)
  
  
  mdatab1=databucket[[1]]
  xb1=databucket[[2]]
  
  xb=c(x[1],xb1,x[m])
  
  mdatab=cbind(mdata[,1],mdatab1,mdata[,m])
  
  
  mdata1b=data.frame(mdatab)
  
  mdatb= data.frame(x=as.numeric(xb), t(mdata1b))
  mdatb= melt(mdatb, id="x")
  
  
  
  rdatabucket=bucket(data[,3:m-1],bb,x[3:m-1],m-2)
  
  rmdatab=cbind(data[,1],rdatabucket[[1]],data[,m])
  rxb=c(x[1],databucket[[2]],x[m])
  
  
  
  
  rmdata1b=data.frame(rmdatab)
  
  rmdatb= data.frame(x=as.numeric(rxb), t(rmdata1b))
  rmdatb= melt(rmdatb, id="x")
  
  
  #miss forest
  
  
  Impmissforest1=missForest(mdatab,maxiter = 10, ntree = 200)   #maxiter = 3, ntree = 20
  Impmissforest=Impmissforest1$ximp
  
  
  
  ymfdat= data.frame(x=as.numeric(xb), t(Impmissforest))
  ymfdat= melt(ymfdat, id="x")
  
  
  #LLR
  
  
  yImpmissforestl1=missForest1(mdatab)   #maxiter = 10, ntree = 100
  yImpmissforestl=yImpmissforestl1$ximp
  
  mfdatl= data.frame(x=as.numeric(xb), t(yImpmissforestl))
  mfdatl= melt(mfdatl, id="x")
  
  
  
  #PACE
  
  
  mlist <- vector("list", nrow(mdatab))
  for (i in 1:nrow(mdatab)) {
    mlist[[i]] <- mdatab[i,]
  }
  
  #mlist[!is.na(mlist)]
  mlisty=lapply(mlist, function(x) x[!is.na(x)])
  
  qq=matrix(rep(xb, each=n), ncol=buck)
  
  mm=mdatab*0
  qw=qq+mm
  
  mlist <- vector("list", nrow(qw))
  for (i in 1:nrow(qw)) {
    mlist[[i]] <- qw[i,]
  }
  
  mlistx=lapply(mlist, function(x) x[!is.na(x)])
  
  
  
  
  res= FPCA(mlisty, mlistx,
            list(dataType='Sparse', error=FALSE,nRegGrid = buck, kernel='epan', verbose=TRUE,userBwCov = 2))
  
  X.pace=fitted.values(res)
  
  pacedat=data.frame(x=xb,t(X.pace))
  pacedat=melt(pacedat,id='x')
  
  
  #   #Imputation with PACE and MISS FOREST
  
  xmis=X.pace
  yImpmissforesttry=pacemissForestb(mdatab)   
  yImpmissforesttry=yImpmissforesttry$ximp
  
  ymfdattry= data.frame(x=as.numeric(xb), t(yImpmissforesttry))
  ymfdattry= melt(ymfdattry, id="x")
  
  
  ####LLR pace
  
  
  yImpmissforesttry1=missForest2b(mdatab)   
  yImpmissforesttry1=yImpmissforesttry1$ximp
  
  ymfdattry1= data.frame(x=as.numeric(xb), t(yImpmissforesttry1))
  ymfdattry1= melt(ymfdattry1, id="x")
  
  
  
  ####### interpolation
  
  
  
  yapprox=function(x,ya,xb){
    yyy=matrix(NA,n,m)
    for(i in 1:n){
      
      yy=ya[i,]
      
      yyy[i,]=c(yy[1],approx(xb,yy,xout=x[3:m-1])$y,yy[buck])
      
    }
    
    return(yyy)
    
  }
  
  ymiss=yapprox(x,Impmissforest,xb)
  
  
  ymfdat= data.frame(x=as.numeric(x), t(ymiss))
  ymfdat= melt(ymfdat, id="x")
  
  
  
  yllr=yapprox(x,yImpmissforestl,xb)
  
  mfdatl= data.frame(x=as.numeric(x), t(yllr))
  mfdatl= melt(mfdatl, id="x")
  
  
  
  ymissp=yapprox(x,yImpmissforesttry,xb)
  
  ymfdattry= data.frame(x=as.numeric(x), t(ymissp))
  ymfdattry= melt(ymfdattry, id="x")
  
  
  
  
  
  yllrp=yapprox(x,yImpmissforesttry1,xb)
  
  ymfdattry1= data.frame(x=as.numeric(x), t(yllrp))
  ymfdattry1= melt(ymfdattry1, id="x")
  
  
  
  
  
  r1=rmse(ymiss,data)
  r2=rmse(yllr,data)
  r3=rmse(ymissp,data)
  r4=rmse(yllrp,data)
  
  
  
  
  ####Spline
  
  
  yspline=function(x,ya,xb){
    
    yyy=matrix(NA,n,m)
    for(i in 1:n){
      
      yy=ya[i,]
      
      yyy[i,]=spline(xb, yy, n =m, method = "fmm",xmin = 0, xmax = 1)$y
      
    }
    
    
    return(yyy)
    
  }
  
  ymiss1=yspline(x,Impmissforest,xb)
  
  
  
  ymfdat1= data.frame(x=as.numeric(x), t(ymiss1))
  ymfdat1= melt(ymfdat1, id="x")
  
  
  
  yllr1=yspline(x,yImpmissforestl,xb)
  
  mfdatl1= data.frame(x=as.numeric(x), t(yllr1))
  mfdatl1= melt(mfdatl1, id="x")
  
  
  
  ymissp1=yspline(x,yImpmissforesttry,xb)
  
  ymfdattry01= data.frame(x=as.numeric(x), t(ymissp1))
  ymfdattry01= melt(ymfdattry01, id="x")
  
  
  
  
  
  yllrp1=yspline(x,yImpmissforesttry1,xb)
  
  ymfdattry11= data.frame(x=as.numeric(x), t(yllrp1))
  ymfdattry11= melt(ymfdattry11, id="x")
  
  
  
  r11=rmse(ymiss1,data)
  r21=rmse(yllr1,data)
  r31=rmse(ymissp1,data)
  r41=rmse(yllrp1,data)
  
  
  
  
  #### model
  
  
  #MODEL
  
  beta=sin(x*2*pi)
  y=data%*%beta/m+rnorm(n)
  ydata=cbind(data,y)
  
  
  truemodel1=pfr(y~fpc(data))
  truemodel1$fitted.values
  t1=rmse(y,data.frame(truemodel1$fitted.values))
  t11=rmse(beta,coef(truemodel1)$value)
  
  
  modelmiss=pfr(y~fpc(ymiss))
  modelmiss$fitted.values
  t2=rmse(y,data.frame(modelmiss$fitted.values))
  t22=rmse(beta,coef(modelmiss)$value)
  
  modelmissp=pfr(y~fpc(ymissp))
  modelmissp$fitted.values
  tp2=rmse(y,data.frame(modelmissp$fitted.values))
  tp22=rmse(beta,coef(modelmissp)$value)
  
  modelmiss1=pfr(y~fpc(ymiss1))
  modelmiss1$fitted.values
  t12=rmse(y,data.frame(modelmiss1$fitted.values))
  t122=rmse(beta,coef(modelmiss1)$value)
  
  
  modelmissp1=pfr(y~fpc(ymissp1))
  modelmissp1$fitted.values
  tp12=rmse(y,data.frame(modelmissp1$fitted.values))
  tp122=rmse(beta,coef(modelmissp1)$value)
  
  
  
  modelllr=pfr(y~fpc(yllr))
  modelllr$fitted.values
  t3=rmse(y,data.frame(modelllr$fitted.values))
  t33=rmse(beta,coef(modelllr)$value)
  
  
  
  modelllrp=pfr(y~fpc(yllrp))
  modelllrp$fitted.values
  tp3=rmse(y,data.frame(modelllrp$fitted.values))
  tp33=rmse(beta,coef(modelllrp)$value)
  
  
  modelllr1=pfr(y~fpc(yllr1))
  modelllr1$fitted.values
  t13=rmse(y,data.frame(modelllr1$fitted.values))
  t133=rmse(beta,coef(modelllr1)$value)
  
  
  
  modelllrp1=pfr(y~fpc(yllrp1))
  modelllrp1$fitted.values
  tp13=rmse(y,data.frame(modelllrp1$fitted.values))
  tp133=rmse(beta,coef(modelllrp1)$value)
  
  
  pmfcoefs = data.frame(grid = x,
                        data= coef(truemodel1)$value,
                        
                        miss_forest= coef(modelmiss)$value,
                        miss_forest_p= coef(modelmissp)$value,
                        LLR= coef(modelllr)$value,
                        LLR_p= coef(modelllrp)$value,
                        
                        smiss_forest= coef(modelmiss1)$value,
                        smiss_forest_p= coef(modelmissp1)$value,
                        sLLR= coef(modelllr1)$value,
                        sLLR_p= coef(modelllrp1)$value,
                        
                        
                        Truth = beta)
  coefs.m = melt(pmfcoefs, id = "grid")
  colnames(coefs.m) = c("grid", "Method", "Value")
  qqmodel1=ggplot(coefs.m, aes(x = grid, y = Value, color = Method, group
                               = Method),width=12,height=6) + geom_path() + theme_bw()
  
  
  mm=matrix(rep(beta,n),n,m, byrow = T)
  mm1=matrix(rep(coef(modelmiss)$value,n),n,m, byrow = T)
  mm2=matrix(rep(coef(modelmissp)$value,n),n,m, byrow = T)
  mm3=matrix(rep(coef(modelllr)$value,n),n,m, byrow = T)
  mm4=matrix(rep(coef(modelllrp)$value,n),n,m, byrow = T)
  
  mms1=matrix(rep(coef(modelmiss1)$value,n),n,m, byrow = T)
  mms2=matrix(rep(coef(modelmissp1)$value,n),n,m, byrow = T)
  mms3=matrix(rep(coef(modelllr1)$value,n),n,m, byrow = T)
  mms4=matrix(rep(coef(modelllrp1)$value,n),n,m, byrow = T)
  
  t000=0
  t111=(rmse(ymiss*mm1,mm*data))
  t222=(rmse(ymissp*mm2,mm*data))
  t333=(rmse(yllr*mm3,mm*data))
  t444=(rmse(yllrp*mm4,mm*data))
  
  ts111=(rmse(ymiss1*mms1,mm*data))
  ts222=(rmse(ymissp1*mms2,mm*data))
  ts333=(rmse(yllr1*mms3,mm*data))
  ts444=(rmse(yllrp1*mms4,mm*data))
  
  
  rall=c(0,r1,r2,r3,r4,r11,r21,r31,r41)
  tall=c(t1,t2,tp2,t12,tp12,t3,tp3,t13,tp13)
  tall1=c(t11,t22,tp22,t122,tp122,t33,tp33,t133,tp133)
  tall2=c(0,t111,t222,t333,t444,ts111,ts222,ts333,ts444)
  tq=c("data","Missforest","Miss forest P","LLR","LLR_P","S Missforest","S Miss forest P","S LLR","S LLR_P")
  tmatrix=matrix(c(tall,tall1,tall2,tq,rall),ncol=5,nrow=9)
  
  return(tmatrix)
  
  
}



ten_interation_nb=function(n,m,o,seeds,buck){
  
  w1=proc.time()
  rr=list(0)
  for(i in 1:length(seeds)){
    
    rr[[i]]=allmodelimputeb(n=n,m=m,o=o,buck=buck,seeds=seeds[i])
    
    
  }
  w2=w1-proc.time()
  
  return(list(rr,w2))
  
}




fit_ten_1000_1=ten_interation_nb(n=1000,m=52,o=.9,seeds = c(1,2,3,4,5,6,7,8,9,10),buck=7)
fit_ten_1000_2=ten_interation_nb(n=1000,m=52,o=.9,seeds = c(1,2,3,4,5,6,7,8,9,10),buck=12)
fit_ten_1000_3=ten_interation_nb(n=1000,m=52,o=.9,seeds = c(1,2,3,4,5,6,7,8,9,10),buck=27)





aa=fit_ten_1000_1
er1=as.numeric(aa[[1]][[1]][,1])+as.numeric(aa[[1]][[2]][,1])+as.numeric(aa[[1]][[3]][,1])+as.numeric(aa[[1]][[4]][,1])+as.numeric(aa[[1]][[5]][,1])+as.numeric(aa[[1]][[6]][,1])+as.numeric(aa[[1]][[7]][,1])+as.numeric(aa[[1]][[8]][,1])+as.numeric(aa[[1]][[9]][,1])+as.numeric(aa[[1]][[10]][,1])
er2=as.numeric(aa[[1]][[1]][,2])+as.numeric(aa[[1]][[2]][,2])+as.numeric(aa[[1]][[3]][,2])+as.numeric(aa[[1]][[4]][,2])+as.numeric(aa[[1]][[5]][,2])+as.numeric(aa[[1]][[6]][,2])+as.numeric(aa[[1]][[7]][,2])+as.numeric(aa[[1]][[8]][,2])+as.numeric(aa[[1]][[9]][,2])+as.numeric(aa[[1]][[10]][,2])
er3=as.numeric(aa[[1]][[1]][,3])+as.numeric(aa[[1]][[2]][,3])+as.numeric(aa[[1]][[3]][,3])+as.numeric(aa[[1]][[4]][,3])+as.numeric(aa[[1]][[5]][,3])+as.numeric(aa[[1]][[6]][,3])+as.numeric(aa[[1]][[7]][,3])+as.numeric(aa[[1]][[8]][,3])+as.numeric(aa[[1]][[9]][,3])+as.numeric(aa[[1]][[10]][,3])
er5=as.numeric(aa[[1]][[1]][,5])+as.numeric(aa[[1]][[2]][,5])+as.numeric(aa[[1]][[3]][,5])+as.numeric(aa[[1]][[4]][,5])+as.numeric(aa[[1]][[5]][,5])+as.numeric(aa[[1]][[6]][,5])+as.numeric(aa[[1]][[7]][,5])+as.numeric(aa[[1]][[8]][,5])+as.numeric(aa[[1]][[9]][,5])+as.numeric(aa[[1]][[10]][,5])
er1=er1/10
er2=er2/10
er3=er3/10
er5=er5/10


er1=matrix(c(er1,er2,er3,er5),9,4,byrow = F)

ert_1000_1=cbind(er1,fit_ten_1000_1[[1]][[1]][,4])





aa=fit_ten_1000_2
er1=as.numeric(aa[[1]][[1]][,1])+as.numeric(aa[[1]][[2]][,1])+as.numeric(aa[[1]][[3]][,1])+as.numeric(aa[[1]][[4]][,1])+as.numeric(aa[[1]][[5]][,1])+as.numeric(aa[[1]][[6]][,1])+as.numeric(aa[[1]][[7]][,1])+as.numeric(aa[[1]][[8]][,1])+as.numeric(aa[[1]][[9]][,1])+as.numeric(aa[[1]][[10]][,1])
er2=as.numeric(aa[[1]][[1]][,2])+as.numeric(aa[[1]][[2]][,2])+as.numeric(aa[[1]][[3]][,2])+as.numeric(aa[[1]][[4]][,2])+as.numeric(aa[[1]][[5]][,2])+as.numeric(aa[[1]][[6]][,2])+as.numeric(aa[[1]][[7]][,2])+as.numeric(aa[[1]][[8]][,2])+as.numeric(aa[[1]][[9]][,2])+as.numeric(aa[[1]][[10]][,2])
er3=as.numeric(aa[[1]][[1]][,3])+as.numeric(aa[[1]][[2]][,3])+as.numeric(aa[[1]][[3]][,3])+as.numeric(aa[[1]][[4]][,3])+as.numeric(aa[[1]][[5]][,3])+as.numeric(aa[[1]][[6]][,3])+as.numeric(aa[[1]][[7]][,3])+as.numeric(aa[[1]][[8]][,3])+as.numeric(aa[[1]][[9]][,3])+as.numeric(aa[[1]][[10]][,3])
er5=as.numeric(aa[[1]][[1]][,5])+as.numeric(aa[[1]][[2]][,5])+as.numeric(aa[[1]][[3]][,5])+as.numeric(aa[[1]][[4]][,5])+as.numeric(aa[[1]][[5]][,5])+as.numeric(aa[[1]][[6]][,5])+as.numeric(aa[[1]][[7]][,5])+as.numeric(aa[[1]][[8]][,5])+as.numeric(aa[[1]][[9]][,5])+as.numeric(aa[[1]][[10]][,5])
er1=er1/10
er2=er2/10
er3=er3/10
er5=er5/10


er1=matrix(c(er1,er2,er3,er5),9,4,byrow = F)

ert_1000_2=cbind(er1,fit_ten_1000_1[[1]][[1]][,4])





aa=fit_ten_1000_3
er1=as.numeric(aa[[1]][[1]][,1])+as.numeric(aa[[1]][[2]][,1])+as.numeric(aa[[1]][[3]][,1])+as.numeric(aa[[1]][[4]][,1])+as.numeric(aa[[1]][[5]][,1])+as.numeric(aa[[1]][[6]][,1])+as.numeric(aa[[1]][[7]][,1])+as.numeric(aa[[1]][[8]][,1])+as.numeric(aa[[1]][[9]][,1])+as.numeric(aa[[1]][[10]][,1])
er2=as.numeric(aa[[1]][[1]][,2])+as.numeric(aa[[1]][[2]][,2])+as.numeric(aa[[1]][[3]][,2])+as.numeric(aa[[1]][[4]][,2])+as.numeric(aa[[1]][[5]][,2])+as.numeric(aa[[1]][[6]][,2])+as.numeric(aa[[1]][[7]][,2])+as.numeric(aa[[1]][[8]][,2])+as.numeric(aa[[1]][[9]][,2])+as.numeric(aa[[1]][[10]][,2])
er3=as.numeric(aa[[1]][[1]][,3])+as.numeric(aa[[1]][[2]][,3])+as.numeric(aa[[1]][[3]][,3])+as.numeric(aa[[1]][[4]][,3])+as.numeric(aa[[1]][[5]][,3])+as.numeric(aa[[1]][[6]][,3])+as.numeric(aa[[1]][[7]][,3])+as.numeric(aa[[1]][[8]][,3])+as.numeric(aa[[1]][[9]][,3])+as.numeric(aa[[1]][[10]][,3])
er5=as.numeric(aa[[1]][[1]][,5])+as.numeric(aa[[1]][[2]][,5])+as.numeric(aa[[1]][[3]][,5])+as.numeric(aa[[1]][[4]][,5])+as.numeric(aa[[1]][[5]][,5])+as.numeric(aa[[1]][[6]][,5])+as.numeric(aa[[1]][[7]][,5])+as.numeric(aa[[1]][[8]][,5])+as.numeric(aa[[1]][[9]][,5])+as.numeric(aa[[1]][[10]][,5])
er1=er1/10
er2=er2/10
er3=er3/10
er5=er5/10


er1=matrix(c(er1,er2,er3,er5),9,4,byrow = F)

ert_1000_3=cbind(er1,fit_ten_1000_1[[1]][[1]][,4])







ert_1000_1
ert_1000_2
ert_1000_3



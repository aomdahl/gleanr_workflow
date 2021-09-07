################################################################################################################################
#Modified 4/26 by Ashton Omdahl
## Copied over from another function since this wasn't working originally.
################################################################################################################################

#helpful code for debugging
if(FALSE)
{
  source("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/quickLoadData.R")
  all <- quickLoadFactorization("Z")
  X <- all$X
  W <- all$W
  option <- all$option
  option$traitSpecificVar <- TRUE
  option$parallel <- FALSE
}

Update_FL <- function(X, W, option, preF = NULL, preL = NULL){
  # number of features - to avoid using T in R
  D      = ncol(X)
  tStart0 = Sys.time()
  
  #Random initialization
  if(!option[['preinitialize']])
  {
    message('Random initialization...')
    FactorM   = matrix(runif(D*(option[['K']] - 1)), nrow = D);
    #Ashton added in- option to include column of 1s
    if(option[['ones_plain']])
    {
      FactorM = cbind(rep(1, D), FactorM);
    } else if(option[['ones_mixed_std']]) {
      message("Using a standard ubiq factor")
      cor_struct <- cor(X)
      #set up directions (-1 or 1) based on most common direction across the board.
      negatives <- apply(cor_struct, 1, function(x) sum((x < 0)))
      positives <- apply(cor_struct, 1, function(x) sum(x > 0) - 1) #-1 because every one is positively correlated with itself.
      best_ref <- positives >= negatives #this provides the most common reference
      #1 accept the most common reference as the directions- just assume this is the dominating effect
      ones <- ifelse(best_ref, 1, -1)
    } else if(option[['ones_mixed_ref']]) {
      cor_struct <- cor(X)
      #set up directions (-1 or 1) based on most common direction across the board.
      negatives <- apply(cor_struct, 1, function(x) sum((x < 0)))
      positives <- apply(cor_struct, 1, function(x) sum(x > 0) - 1) #-1 because every one is positively correlated with itself.
      best_ref <- positives >= negatives #this provides the most common reference
      #2 find some trait that matches this as closely as possible, let that be the reference and use the directions relative to that one.
      cor_ref <- cor_struct >= 0
      reference_trait <- which.max(apply(cor_ref, 2, function(x) sum(x == best_ref)))
      ones <- ifelse(cor_struct[,reference_trait] >= 0, 1, -1)
    }	else if(option[['ones_eigenvect']]) {
      message("1st column based on direction of svd of cor")
      cor_struct <- cor(X)
      svd <- svd(cor_struct, nu = 1) #fortunately its symmetric, so  U and V are the same here!
      ones <- sign(svd$u)
    } else {
      #FactorM   = matrix(runif(D*(option[['K']])), nrow = D);
      ones = matrix(runif(D*(option[['K']])), nrow = D)[,1]
    }
    FactorM = cbind(ones, FactorM);
  } else #you pre-provide the first F.
  {
    FactorM   = preF;
    message("initialized with PCA!")
  }

  # First round of optimization
  message('Start optimization ...')
  message(paste0('K = ', (option[['K']]), '; alpha1 = ', (option[['alpha1']]),'; lambda1 = ', (option[['lambda1']])));
  
  #Need to adjust for first run- no reweighting (obvi.)
  og_option <- option[['reweighted']]
  option[['reweighted']] <- FALSE
  
  if(option[['parallel']]) #This is not working at all. Can't tell you why. But its not. Need to spend some time debugging at some point.
  {
    L = fit_L_parallel(X, W, FactorM, option, formerL = preL); #preL is by default Null, unless yo specify!
  } else
  {
    L = fit_L(X, W, FactorM, option, formerL = preL);
  }
  
  option[['reweighted']] <- og_option
  objective = c(NA, compute_obj(X, W, L, FactorM, option));
  objective_change = c(1, 1);
  old_change = 1;
  F_old = FactorM; 
  
  # if L is empty, stop
  non_empty_l = which(apply(L, 2, function(x) sum(x!=0)) > 0)
  if(length(non_empty_l) == 0){
    message('Finished');
    message('L is completely empty, alpha1 too large')
    FactorM = NULL;
    F_sparsity = 1;
    L_sparsity = 1;
    factor_corr = 1;
    Nfactor = 0;
    return()
  }
  

    trait.var <- matrix(NA, option[['iter']], ncol(X))

  for (iii in seq(1, option[['iter']])){
    
    ## update F
    if(iii == 1)
    {
      og_option <- option[['reweighted']]
      option[['reweighted']] <- FALSE
      FactorM = fit_F(X, W, L, option, FactorM);
      option[['reweighted']] <- og_option
    } else
    {
      FactorM = fit_F(X, W, L, option, FactorM); #by iter 3 really slows down, due to the L1 requirements. Yea this won't do....
    }
    message(paste0("Currently on iteration ", iii))
    #get the factor specific variance....
    if(option$traitSpecificVar)
    {
      trait.var[iii,] <- FactorM$r.v
      FactorM <- FactorM$mat
    }
    # if number of factors decrease because of empty factor, the change in ||F||_F = 100
    if(ncol(FactorM) != ncol(F_old)){
      F_change = 100
    }else{
      F_change = norm(FactorM - F_old, 'F') / ncol(FactorM)
    }
    F_old = FactorM;
    
    non_empty_f = which(apply(FactorM, 2, function(x) sum(x!=0)) > 0)
    if(length(non_empty_f) == 0){
      message('Finished');
      message('F is completely empty, lambda1 too large')
      F_sparsity = 1;
      L_sparsity = 1;
      factor_corr = 1;
      Nfactor = 0;
      return();
    }
    colnames(FactorM) = seq(1, ncol(FactorM));
    
    ## update L
    if(option[['parallel']])
    {
      L  = fit_L_parallel(X, W, FactorM, option, L);
    }else{
      L  = fit_L(X, W, FactorM, option, L, r.v = trait.var[iii,]); #the l1 one requirement is making things tough here.....
    }

    # if L is empty, stop
    non_empty_l = which(apply(L, 2, function(x) sum(x!=0)) > 0)
    if(length(non_empty_l) == 0){
      message('Finished');
      message('L is completely empty, alpha1 too large')
      FactorM = NULL;
      F_sparsity = 1;
      L_sparsity = 1;
      factor_corr = 1;
      Nfactor = 0;
      return()
    }
    colnames(L) = seq(1, ncol(L));
    
    # align the two matrices	
    non_empty = intersect(non_empty_l,non_empty_f);
    L = as.matrix(as.data.frame(L[, non_empty]));
    FactorM  = as.matrix(as.data.frame(FactorM[, non_empty]));
    
    
    # collect sparsity in L and F
    L_sparsity = sum(abs(L) < 1e-5) / ncol(L) / nrow(L);
    F_sparsity = sum(abs(FactorM) < 1e-5) / ncol(FactorM) / nrow(FactorM);
    Nfactor = length(non_empty);
    
    # change in the objective -- should always be negative
    obj_updated = compute_obj(X, W, L, FactorM, option);
    obj_change = obj_updated - objective[length(objective)];
    objective = c(objective, obj_updated);
    print(objective)
    objective_change = c(objective_change, obj_change);
    
    if(option[['disp']]){
      cat('\n')
      message(paste0('Iter', iii, ':'))
      message(paste0('Objective change = ', obj_change))
      message(paste0('Frobenius norm of (updated factor matrix - previous factor matrix) / number of factors  = ', F_change));
      message(paste0('Loading Sparsity = ', L_sparsity, '; Factor sparsity = ', F_sparsity, '; ', Nfactor, ' factors remain')); 
      cat('\n')
    }
    
    
    # converge if: 1). Change of the values in factor matrix is small. ie. The factor matrix is stable. 2). Change in the objective function becomes small; 3). reached maximum number of iterations
    if(option[['convF']] >= F_change){
      message(paste0('Factor matrix converges at iteration ', iii));
      break
    }
    
    oc1 = abs(objective_change[length(objective_change)])
    if(option[['convO']] >= oc1){
      message(paste0('Objective function converges at iteration ', iii));
      break
    }
    if(iii == option[['iter']]){
      message('Reached maximum iteration.');
      break
    }
    
  }
  
  tEnd0 = Sys.time()
  cat('\n')
  message('Total time used for optimization: ');
  message(tEnd0 - tStart0);
  cat('\n')
  
  # return L, F, sparsity in L and F, number of factors -- could be different from K!
  return(list(FactorM, L, L_sparsity, F_sparsity, Nfactor, objective, trait.var))
}

#9/07 testing notes
#tried regular, works
#trying now with factor-specific variance, shall see.

#ran
#t2 <- Update_FL(X,W, option)
library(ggplot)
#res <- data.frame(t2[[7]]) %>% mutate("iter" = c(1,2,3,4)) %>% pivot_longer(cols = paste0("X", 1:55))
#ggplot(res, aes(x = iter, y  = value, color = name)) + geom_line()
option$alpha1 <- 0
option$lambda1 <- 0
option$iter <- 10
t3 <- Update_FL(X,W, option)

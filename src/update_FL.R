################################################################################################################################
#Modified 4/26 by Ashton Omdahl
## Copied over from another function since this wasn't working originally.
################################################################################################################################

#helpful code for debugging
if(FALSE)
{
  #source("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/quickLoadData.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/quickLoadData.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_F.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_L.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/compute_obj.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/plot_functions.R")
  all <- quickLoadFactorization("Z", "MARCC")
  X <- all$X
  W <- all$W
  option <- all$option
  option$traitSpecificVar <- TRUE
  option$parallel <- FALSE
  
  #subset for faster running....
  X <- X[1:1000, 1:10]
  W <- W[1:1000, 1:10]
  option$K <- 5
}

      #fast matrix correlation help
library(coop)
cor2 <- function(x) {
1/(NROW(x)-1) * crossprod(scale(x, TRUE, TRUE))
}
      

Update_FL <- function(X, W, option, preF = NULL, preL = NULL){
  # number of features - to avoid using T in R
  D = ncol(X)
  tStart0 = Sys.time()
  FactorM = NULL 
  #Random initialization of F
  if(option$l_init != "")
  {
      message("Initializing L rather than F.")
      nsnps = nrow(X)
      library(irlba)

      L   = matrix(runif(nsnps*(option[['K']] - 1)), nrow = nsnps);
      if(option$l_init == 'ones_plain')
      {
        message("Initializing ubiquitous factor with all ones...")
        FactorM = cbind(rep(1, nsnps), FactorM);

      }	else if(option$l_init == 'ones_eigenvect') {
        message("1st column based on direction of svd of cor")
        cor_struct <- cor2(t(X))
        svd <- irlba(cor_struct, 1) #fortunately its symmetric, so  U and V are the same here!
        ones <- sign(svd$u)
      } else if(option$l_init == 'plieotropy')
      {
        message("1st column based svd of cor(|Z|), since plieotropy has no direction.")
        cor_struct <- cor2(abs(t(X)))
        svd <- irlba(cor_struct, 1) #fortunately its symmetric, so  U and V are the same here!
        ones <- svd$u
      } else {
        #FactorM   = matrix(runif(D*(option[['K']])), nrow = D);
        ones = matrix(runif(nsnps*(option[['K']])), nrow = D)[,1]
      }
      L = cbind(ones, L);
      # First round of optimization
      message('Start optimization ...')
      message(paste0('K = ', (option[['K']]), '; alpha1 = ', (option[['alpha1']]),'; lambda1 = ', (option[['lambda1']])));
      
      FactorM = fit_F(X, W, L, option)
      
  } else{
    if(!option[['preinitialize']])
    {
      FactorM   = matrix(runif(D*(option[['K']] - 1)), nrow = D);
      #Ashton added in- option to include column of 1s
      if(option[['f_init']] == 'ones_plain')
      {
        message("Initializing ubiquitous factor with all ones...")
        FactorM = cbind(rep(1, D), FactorM);
      }	else if(option[['f_init']] == 'ones_eigenvect') {
        message("1st column based on direction of svd of cor")
        cor_struct <- cor(X)
        svd <- svd(cor_struct, nu = 1) #fortunately its symmetric, so  U and V are the same here!
        ones <- sign(svd$u)
      } else if(option[['f_init']] == 'plieotropy')
      {
        message("1st column based svd of cor(|Z|), since plieotropy has no direction.")
        cor_struct <- cor(abs(X))
        svd <- svd(cor_struct, nu = 1) #fortunately its symmetric, so  U and V are the same here!
        ones <- svd$u
      } else {
        #FactorM   = matrix(runif(D*(option[['K']])), nrow = D);
        ones = matrix(runif(D*(option[['K']])), nrow = D)[,1]
      }
      FactorM = cbind(ones, FactorM);
    } else #you pre-provide the first F.
    {
      FactorM   = preF;
    }
    if(option$posF)
    {
      message("Performing semi-non-negative factorization today....")
      FactorM = abs(FactorM)
    }
  }
  
# First round of optimization
message('Start optimization ...')
message(paste0('K = ', (option[['K']]), '; alpha1 = ', (option[['alpha1']]),'; lambda1 = ', (option[['lambda1']])));

#Need to adjust for first run- no reweighting (obvi.)
og_option <- option[['reweighted']]
option[['reweighted']] <- FALSE
if(option[['ncores']] > 1) #This is not working at all. Can't tell you why. But its not. Need to spend some time debugging at some point.
{
  print("Fitting L in parallel")	
  L = fit_L_parallel(X, W, FactorM, option, formerL = preL); #preL is by default Null, unless yo specify!
}
else
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
  #Alternatively, if L has only the one factor and we are doing the unwighted ubiq mode...
  if(option[["fixed_ubiq"]] & (ncol(L) == 1  | ncol(FactorM) == 1))
  {
    message("Only the ubiquitous factor remains, when we are removing any L1 prior on it")
    FactorM = NULL;
    F_sparsity = 1;
    L_sparsity = 1;
    factor_corr = 1;
    Nfactor = 0;
    return()
  }
  
  
  trait.var <- matrix(NA, option[['iter']], ncol(X))
  lambda_track <- c(option[['lambda1']])
  alpha_track <- c(option[['alpha1']])
  for (iii in seq(1, option[['iter']])){
    message(paste0("Currently on iteration ", iii))
    ## update F
    if(option$autofit > 0 &  iii > 1)
    {
        #F
        message(paste0("Current lambda: ", option$lambda1))
        option[['lambda1']] <- MAPfitLambda(FactorM,option$autofit, option)
        message(paste0("Updated lambda: ", option[['lambda1']]))
        #L
        message(paste0  ("Current alpha: ", option$alpha1))
        option[['alpha1']] <- MAPfitAlpha(L,option$autofit,option)
        message(paste0("Updated alpha: ", option[['alpha1']]))
        lambda_track <- c(lambda_track, option[['lambda1']])
        alpha_track <- c(alpha_track, option[['alpha1']])
    }

    if(iii == 1)
    {
      og_option <- option[['reweighted']]
      option[['reweighted']] <- FALSE
      FactorM = fit_F(X, W, L, option, FactorM);
      option[['reweighted']] <- og_option
    } else {
      FactorM = fit_F(X, W, L, option, FactorM); #by iter 3 really slows down, due to the L1 requirements. Yea this won't do....
    }
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
    if(length(non_empty_f) == 0 | (non_empty_f[1] == 1 & option$fixed_ubiq & (length(non_empty_f) == 1))){
      message('Finished');
      message('F is completely empty or loaded only on ubiquitous factor, lambda1 too large')
      F_sparsity = 1;
      L_sparsity = 1;
      factor_corr = 1;
      Nfactor = 0;
      if(length(non_empty_f) == 0) {
        return()
        }
      else {
        L = matrix(0, ncol(FactorM), 1)
        break
        }
    }
    colnames(FactorM) = seq(1, ncol(FactorM));
    
    ## update L
    #message("Fitting L....")
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
    print(paste0("Objective change: ", oc1))
    print(paste0("Target objective change: ",option[['convO']]))
    if(oc1 < option[['convO']]){
      message("Objective function change threshold achieved!")
      message(paste0('Objective function converges at iteration ', iii));
      break
    }
    if(iii == option[['iter']]){
      message('Reached maximum iteration.');
      break
    }
  }
  retlist <- list("F" = FactorM, "L" = L, "L_sparsity" = L_sparsity, "F_sparsity" = F_sparsity, "K" = Nfactor, "obj" = objective, "study_var" = trait.var)
  tEnd0 = Sys.time()
  cat('\n')
  message('Total time used for optimization: ');
  message(tEnd0 - tStart0);
  cat('\n')
  # return L, F, sparsity in L and F, number of factors -- could be different from K!
    return(retlist) #F, L, l sparsity ,f_sparsity, K, obj, study_var, 
}


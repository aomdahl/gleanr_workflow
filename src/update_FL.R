################################################################################################################################
#Modified 4/26 by Ashton Omdahl
## Copied over from another function since this wasn't working originally.
################################################################################################################################

#helpful code for debugging
if(FALSE)
{
  #source("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/quickLoadData.R")
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/quickLoadData.R")
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_F.R")
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_L.R")
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/compute_obj.R")
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/plot_functions.R")
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
ZERO_THRESH = 1e-5

  #fast matrix correlation help
library(coop)
cor2 <- function(x) {
1/(NROW(x)-1) * crossprod(scale(x, TRUE, TRUE))
}
      
#This function uses OLS to get a good estimate of what the maximum sparsity space is across all parameters.
#@param burn.in- how many iterations to go
#@return V_burn: the burn-in estimate of V
#@return U_burn: the burn-in estimate of U
#@return max_sparsity_params: a list indexed by iteration containing the list of max alphas and lambdas (i.e. list[[i]]$alpha))
defineSparsitySpace <- function(X, W, option, burn.in = 5)
{
  #Perform burn.in # of iterations with no sparsity (OLS) and calculate the maximum sparsity value at each step.
  new.options <- option
  new.options$burn.in <- burn.in
  new.options$regression_method <- "OLS"
  new.options$calibrate_sparsity <- TRUE
  new.options$actively_calibrating_sparsity <- TRUE
  #Record those sparsities PARAMETERS for each iteration for both L and F, and then return for downstream analysis.
  param.space <- list()
  #for testing purposes
  Xint <- as.matrix(X)
  Wint <-as.matrix(W)
  V <- initV(Xint,Wint, new.options)
  for(i in 1:burn.in)
  {
    U.dat <- FitUWrapper(Xint,Wint,V, new.options)
    V.dat <- FitVWrapper(Xint, Wint, U.dat$U, new.options, FactorM);
    param.space[[i]] <- list("alpha" = U.dat$sparsity_space, "lambda"=V.dat$sparsity_space)
    V <- V.dat$V
  }
  #Upon return, we want the max across all as the maximums.
  #Maybe set the min and max of the distribution as the upper and lower bounds, and then pick some points in between?
  #That would guarantee at least one row or one column would be entirely empty.
  return(list("V_burn" = V, "U_burn"=U.dat$U, "max_sparsity_params"=param.space))
}

#Function to initialize V as desired
initV <- function(X,W,option, preV = NULL)
{
  D = ncol(X)
  if(!option[['preinitialize']])
  {
    V = matrix(runif(D*(option[['K']] - 1)), nrow = D);

    if(option[['f_init']] == 'ones_plain')
    {
      message("Initializing ubiquitous factor with all ones...")
      V = cbind(rep(1, D), V);
    }	else if(option[['f_init']] == 'ones_eigenvect') {
      message("1st column based on direction of svd of cor")
      cor_struct <- cor2(X)
      svd <- svd(cor_struct, nu = option$K) #fortunately its symmetric, so  U and V are the same here!
      ones <- sign(svd$u[,1])
      V <- svd$u[,2:option$K]
    } else if(option[['f_init']] == 'plieotropy')
    {
      message("1st column based svd of cor(|Z|), since plieotropy has no direction.")
      cor_struct <- cor(abs(X))
      svd <- svd(cor_struct, nu = 1) #fortunately its symmetric, so  U and V are the same here!
      ones <- svd$u
    } else {
      ones = matrix(runif(D*(option[['K']])), nrow = D)[,1]
    }
    V = cbind(ones, V);
  } else #you pre-provide the first F.
  {
    V   = preV;
  }
  if(option$posF)
  {
    message("Performing semi-non-negative factorization today....")
    V = abs(V)
  }
  V
}

#Function to inialize U if specified
initU <- function(X,W,option, prevU = NULL)
{
  message("Initializing L rather than F.")
  nsnps = nrow(X)
  library(irlba)
  U   = matrix(runif(nsnps*(option[['K']] - 1)), nrow = nsnps);
  if(option$l_init == 'ones_plain')
  {
    message("Initializing ubiquitous factor with all ones...")
    ones = rep(1, nsnps)
    
  }	else if(option$l_init == 'ones_eigenvect') {
    message("1st column based on direction of svd of cor")
    cor_struct <- cor2(t(X))
    svd <- irlba(cor_struct, 1) #fortunately its symmetric, so  U and V are the same here!
    ones <- sign(svd$u)
  } else if(option$l_init == 'plieotropy') {
    message("1st column based svd of cor(|Z|), since plieotropy has no direction.")
    cor_struct <- cor2(abs(t(X)))
    svd <- irlba(cor_struct, 1) #fortunately its symmetric, so  U and V are the same here!
    ones <- svd$u
  } else {
    ones = matrix(runif(nsnps*(option[['K']])), nrow = D)[,1]
  }
  cbind(ones, L);
  
}

#helper code to clean things up, just keep track of those relevant metrics,
UpdateTrackingParams <- function(storage.unit, X,W,U,V,option, sparsity.thresh = 1e-5,lambda_step = NULL,alpha_step = NULL)
{
  if(is.null(storage.unit))
  {
    objective = c(NA, compute_obj(X, W, L, FactorM, option));
    objective_change = c(NA);
    sto.obj <-  list("U" = U, "V" = V, "K" = ncol(U), "obj" = objective, "obj_change" = objective_change,
                     "V_sparsities" = c(), "U_sparsities" = c(), "autofit_lambda" = c(), "autofit_alpha"=c(), "mse" = c())
  }
    # collect sparsity in L and F
  
  sto.obj$U_sparsities = c(sto.obj$U_sparsities, sum(abs(U) < sparsity.thresh) / (ncol(U) * nrow(U)));
  sto.obj$V_sparsities = c(sto.obj$V_sparsities, sum(abs(V) < sparsity.thresh) / (ncol(V) * nrow(V)));
  sto.obj$K = ncol(U);
    
    # change in the objective -- should always be negative
    obj_updated = compute_obj(X, W, U, V, option);
    obj_change = obj_updated - storage.unit$obj[length(storage.unit$obj)];
    storage.unit$obj = c(storage.unit$obj, obj_updated);
    storage.unit$obj_change = c(storage.unit$obj_change, obj_change);
    #Update the fit:
    sto.obj$mse <- c(sto.obj$mse, calcMSE(X,W,V, U))
    
  #Update the sparsity paramters, regardless of if they change or not.
    if(is.null(alpha_step) & is.null(lambda_step))
    {
      lambda_step <- option$lambda1
      alpha_step <- option$alpha1
    }
      sto.obj$autofit_lambda <- c(sto.obj$autofit_lambda,lambda_step)
      sto.obj$autofit_alpha <- c(sto.obj$autofit_alpha,alpha_step)
  sto.obj
}

#This implements the update step for MAP fits, and includes settings if we arejust fixing one and learning the other (both options)
UpdateSparsityMAPAutofit <- function(iter, U, V, option)
{
  opt.ret <- option
  lambda.prev <- opt.ret[['lambda1']] 
  if(opt.ret$swap) {
    lambda.prev <- opt.ret[['alpha1']] 
  }
  if(opt.ret$autofit != -1 &  iii > 1)
  {
    #F
    message(paste0("Current lambda: ", opt.ret$lambda1))
    opt.ret[['lambda1']] <- MAPfitLambda(V,opt.ret$autofit, opt.ret)
    message(paste0("Updated lambda: ", opt.ret[['lambda1']]))
    #L
    message(paste0  ("Current alpha: ", opt.ret$alpha1))
    opt.ret[['alpha1']] <- MAPfitAlpha(U,opt.ret$autofit,opt.ret)
    message(paste0("Updated alpha: ", opt.ret[['alpha1']]))
    
  }
  if(opt.ret$autofit != -1 & !is.na(opt.ret$fix.alt.setting)) #opt.ret to fix one and modulate the other, given we are doing autofit.
  {
    if(iter < (opt.ret$fix.alt.setting * opt.ret$iter))
    {
      if(opt.ret$swap)
      {
        message("Lambda setting fixed at 0")
        opt.ret[['lambda1']] <- 1e-10
        
      }else
      {
        message("Alpha setting fixed at 0")
        opt.ret[['alpha1']] <- 1e-10
      }
      
    }else
    {
      if(opt.ret$swap)
      {
        message("Alpha setting fixed now")
        print(lambda.prev)
        opt.ret[['alpha1']] <- lambda.prev
        
      }else
      {
        message("Lambda setting fixed now")
        print(lambda.prev)
        opt.ret[['lambda1']] <- lambda.prev
      }
      
    }
  }
  opt.ret
}


#Checking results along the way
CheckUEmpty <- function(U) 
{
  non_empty_u = which(apply(U, 2, function(x) sum(x!=0)) > 0) #Truly all 0
  if(length(non_empty_l) == 0){
    message('Finished');
    message('L is completely empty, alpha1 too large')
    #FactorM = NULL;
    #F_sparsity = 1;
    #L_sparsity = 1;
    #factor_corr = 1;
    #Nfactor = 0;
    return(TRUE)
  }
  FALSE
}

CheckVEmpty <- function(V)
{
  non_empty_v = which(apply(V, 2, function(x) sum(x!=0)) > 0)
  #CHange: only ending if all the factors are empty. This is basically impossible. We proceed if the first factor is still valid...
  if(length(non_empty_v) == 0){  #| (non_empty_v[1] == 1 & option$fixed_ubiq & (length(non_empty_f) == 1))){
    updateLog('Finished', option);
    updateLog('F is completely empty or loaded only on ubiquitous factor, lambda1 too large', option)
    return(TRUE)
  } 
  return(FALSE)
}


AlignFactorMatrices <- function(U, V)
{
  non_empty_v = which(apply(V, 2, function(x) sum(x!=0)) > 0)
  non_empty_u = which(apply(U, 2, function(x) sum(x!=0)) > 0)
  non_empty = intersect(non_empty_u,non_empty_v);
  U = as.matrix(as.data.frame(U[, non_empty]));
  V  = as.matrix(as.data.frame(V[, non_empty]));
  return(list("U"=U, "V"=V))
}
# converge if: 1). Change of the values in factor matrix is small. ie. The factor matrix is stable. 2). Change in the objective function becomes small; 3). reached maximum number of iterations

ConvergenceConditionsMet <- function(U,V,tracker,option)
{
  #TODO double check tracker$V tracks with the old
  #Condition 1: Change in V is small
  V_change = norm(V - tracker$V, 'F') / ncol(V)
  if(option[['convF']] >= V_change){
    updateLog(paste0('Factor matrix converges at iteration ', iii), option);
    return(TRUE)
  }
  #option 2: objective change is small
  oc1 = abs(tracker$objective_change[length(tracker$objective_change)])
  updateLog(paste0("Objective change: ", oc1), option)
  if(oc1 < as.numeric(option[['convO']])){
    updateLog(("Objective function change threshold achieved!"), option)
    updateLog(paste0('Objective function converges at iteration ', iii), option);
    return(TRUE)
  }
  #Option 3- maxi number of iter.
  if(iii == option[['iter']]){
    message('Reached maximum iteration.');
    return(TRUE)
  }
  FALSE
}

#Main workhorse function
Update_FL <- function(X, W, option, preV = NULL, preL = NULL){
  # number of features - to avoid using T in R
  D = ncol(X)
  tStart0 = Sys.time()
  FactorM = NULL 
  #Random initialization of F
  if(option$l_init != "")
  {
    U <- initU(X,W,option,preU=preL)
    V = fit_F(X, W, L, option)$V
      
  } else{ #initialize by F as normal.
    
    V <- initV(X,W,option,preV=preF)
  }  
  
  message("")
  message('Start optimization ...')
  message(paste0('K = ', 
                 (option[['K']]), '; alpha1 = ', round(option[['alpha1']], digits =  4),
                 '; lambda1 = ', round(option[['lambda1']], digits = 4)))

  #Case where you use previous iteration estimates to inform next iteration..
  og_option <- option[['carry_coeffs']]
  option[['carry_coeffs']] <- FALSE
  L = FitUWrapper(X,W,FactorM, option)
  L <- L$U
  option[['carry_coeffs']] <- og_option
  
  #Start tracking stats
  tracking.data <- UpdateTrackingParams(NULL, X,W,L,FactorM,option)
  F_old = tracking.data$V; 
  
  #If U is already empty, than we know that the matrix is too sparse, and we should just end there.
  if(CheckUEmpty(L)) {return(tracking.data)}

  #If we are doing a measure of per-trait variance over time...
  trait.var <- matrix(NA, option[['iter']], ncol(X))
 
  for (iii in seq(1, option[['iter']])){ #start the iterations (although technically, we have done 1 already.)
    message(paste0("Currently on iteration ", iii))
    
    ## If we are doing an autofit setting....
    if(option$autofit) {
      option <- UpdateSparsityMAPAutofit(iii, L,FactorM, option)
    }
    FactorM = fit_F(X, W, L, option, FactorM); #by iter 3 really slows down, due to the L1 requirements. Yea this won't do....
    
    #get the factor specific variance....
    if(option$traitSpecificVar)
    {
      trait.var[iii,] <- FactorM$r.v
      FactorM <- FactorM$mat
    }
    # if number of factors decrease because of empty factor, the change in ||F||_F = 100
    
    #Tracking change in F....
    if(CheckVEmpty(FactorM)) {return(UpdateTrackingParams(tracking.data, X,W,L,FactorM,option))}
    colnames(FactorM) = seq(1, ncol(FactorM));
    
    ## update L
    L <- FitUWrapper(X,W,FactorM, option,r.v = trait.var[iii,])$U

    # if L is empty, stop
    if(CheckUEmpty(L)) {return(UpdateTrackingParams(tracking.data, X,W,L,FactorM,option))} #Update and end.
    colnames(L) = seq(1, ncol(L));
    
    # align the two matrices	
    updated.mats <- AlignFactorMatrices(L, FactorM); L <- updated.mats$U; FactorM <- updated.mats$V

   
    
    
    if(option$V > 0){
      cat('\n')
      updateLog(paste0('Iter', iii, ':'), option)
      updateLog(paste0('Objective change = ', obj_change), option)
      updateLog(paste0('Frobenius norm of (updated factor matrix - previous factor matrix) / number of factors  = ', F_change), option);
      updateLog(paste0('Loading Sparsity = ', L_sparsity, '; Factor sparsity = ', F_sparsity, '; ', Nfactor, ' factors remain'), option); 
      cat('\n')
    }
    
    #pickup here...
    if(ConvergenceConditionsMet(L, FactorM, tracking.data, option))
    {
      cat('\n')
      updateLog(paste0('Total time used for optimization: ',  round(difftime(Sys.time(), tStart0, units = "mins"), digits = 3), ' min'), option);
      cat('\n')
     return(UpdateTrackingParams(tracking.data, X,W,L,FactorM,
                                            option, lambda_step = lambda_track,alpha_step = alpha_track))
    }else{
      tracking.data <- UpdateTrackingParams(tracking.data, X,W,L,FactorM,
                                            option, lambda_step = lambda_track,alpha_step = alpha_track)
    }
   
    
  }
}


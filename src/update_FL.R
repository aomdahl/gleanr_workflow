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
#ZERO_THRESH = 1e-5
#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/sparsity_scaler.R")
#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/pve.R")
  #fast matrix correlation help

cor2 <- function(x) {
  library(coop)
1/(NROW(x)-1) * crossprod(scale(x, TRUE, TRUE))
}
      
#This function uses OLS to get a good estimate of what the maximum sparsity space is across all parameters.
#@param burn.in- how many iterations to go
#@return V_burn: the burn-in estimate of V
#@return U_burn: the burn-in estimate of U
#@return max_sparsity_params: a list indexed by iteration containing the list of max alphas and lambdas (i.e. list[[i]]$alpha))
DefineSparsitySpaceInit <- function(X, W, option, burn.in = 5)
{
  #Perform burn.in # of iterations with no sparsity (OLS) and calculate the maximum sparsity value at each step.
  new.options <- option
  new.options$burn.in <- burn.in
  new.options$regression_method <- "OLS"
  #Record those sparsities PARAMETERS for each iteration for both L and F, and then return for downstream analysis.
  param.space <- list()
  #for testing purposes
  Xint <- as.matrix(X)
  Wint <-as.matrix(W)
  V <- initV(Xint,Wint, new.options)
  for(i in 1:burn.in)
  {
    print(i)
    #U.dat <- FitUWrapper(Xint,Wint,V, new.options)
    U.dat <- DefineSparsitySpace(Xint, Wint, V, "U", new.options, fit = "OLS")
    new.options$K = ncol(U.dat$U) #Update if it has changed
    #If we have ones with NA, they need to get dropped
    #if we have columns with NAs here, we want them gone now so it doesn't jank up downstream stuff.
    #V.dat <- FitVWrapper(Xint, Wint, U.dat$U, new.options, V);
    V.dat <- DefineSparsitySpace(Xint, Wint,U.dat$U, "V", new.options, fit = "OLS")
    if(i > 3)
    {
      #this might be high but oh well...
      #HERE
      #drops <- CheckLowPVE(X,W,V.dat$V) #redundant with other code
      drops <- c()
      print(drops)
      if(length(drops) > 0)
      {
        n <- DropSpecificColumns(drops, U.dat$U, V.dat$V)
        U.dat$U <- n$U
        new.options$K <- ncol(U.dat$U)
        V.dat <- DefineSparsitySpace(Xint, Wint,U.dat$U, "V", new.options, fit = "OLS")
      }
      
    }
    param.space[[i]] <- list("alpha" = U.dat$sparsity_space, "lambda"=V.dat$sparsity_space)
    V <- V.dat$V
  }
  
  #Upon return, we want the max across all as the maximums.
  #Maybe set the min and max of the distribution as the upper and lower bounds, and then pick some points in between?
  #That would guarantee at least one row or one column would be entirely empty.
  return(list("V_burn" = V, "U_burn"=U.dat$U, "max_sparsity_params"=param.space, "new.k" =new.options$K))
}

#TODO: get the organization right- this is the internal function that does the work, Init version is just a wrapper around it.
#    U.dat <- DefineSparsitySpace(Xint, Wint, V, "U", new.options, fit = "OLS")
#Get the sparsity parameters associated with the regression step to learn "loading"
DefineSparsitySpace <- function(X,W,fixed,loading, option, fit = "None")
{
  new.options <- option
  new.options$regression_method <- fit
  new.options$actively_calibrating_sparsity <- TRUE
  if(loading == "V")
  {
    free.dat <- FitVWrapper(X, W, fixed, new.options);
    #HERE
  }else
  {
    red.cols <- c(1)
    while(!is.null(red.cols))
    {
      free.dat <- FitUWrapper(X,W, fixed, new.options)
      red.cols <- free.dat$redundant_cols
      if(!is.null(red.cols)){
        fixed <- fixed[,-red.cols]
        new.options$K <- ncol(fixed)
      }
      
    }
  }
  if(fit == "None")
  {
    return(free.dat$sparsity_space)
  }else{
    return(free.dat)
  }
  
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
    message("Using an initialized V")
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
UpdateTrackingParams <- function(sto.obj, X,W,U,V,option, sparsity.thresh = 1e-5)
{
  if(is.null(sto.obj))
  {
    sto.obj <-  list("V" = V, "U" = U, "initK" = option$K, "K" = ncol(U), "obj" = c(NA), "obj_change" = c(),
                     "V_sparsities" = c(), "U_sparsities" = c(), "autofit_lambda" = c(), "autofit_alpha"=c(), "mse" = c(),
                     "V_change" = c(), "U_change" = c())
  }
    # collect sparsity in L and F
  
  #sto.obj$U_sparsities = c(sto.obj$U_sparsities, sum(abs(U) < sparsity.thresh) / (ncol(U) * nrow(U)));
  #sto.obj$V_sparsities = c(sto.obj$V_sparsities, sum(abs(V) < sparsity.thresh) / (ncol(V) * nrow(V)));
  sto.obj$U_sparsities = c(sto.obj$U_sparsities, matrixSparsity(U, option$K)); #The sparsity count should be with respect to the initial value..
  #TODO: fix the relative sparsity count thuing...
  sto.obj$V_sparsities = c(sto.obj$V_sparsities,matrixSparsity(V, option$K))
  sto.obj$K = ncol(U);
  sto.obj$mse <- c(sto.obj$mse, norm(W*X-W*(U%*%t(V)), type = "F")/(nrow(X) * ncol(X)))
    # change in the objective -- should always be negative
    obj_updated = compute_obj(X, W, U, V, option);
    obj_change = sto.obj$obj[length(sto.obj$obj)] - obj_updated;
  if(is.na(sto.obj$obj[1]))
  {
    sto.obj$obj[1] <- obj_updated
  }
  sto.obj$obj = c(sto.obj$obj, obj_updated);
  sto.obj$obj_change = c(sto.obj$obj_change, obj_change);
    #Update the fit:

  sto.obj$V_change = c(sto.obj$V_change, MatrixChange(V, sto.obj$V))
  sto.obj$U_change = c(sto.obj$U_change, MatrixChange(U, sto.obj$U))
  #Update the sparsity paramters, regardless of if they change or not.
      sto.obj$autofit_lambda <- c(sto.obj$autofit_lambda,option$lambda1)
      sto.obj$autofit_alpha <- c(sto.obj$autofit_alpha,option$alpha1)
      
      #LAST STEP- update the new U and V
      sto.obj$U <- U
      sto.obj$V <- V
      sto.obj$PVE <-PercentVarEx(as.matrix(X)*as.matrix(W), v = V)
  sto.obj
}

#This implments the update step for MAP fits, and includes settings if we arejust fixing one and learning the other (both options)
UpdateSparsityMAPAutofit <- function(iter, U, V, option)
{
  opt.ret <- option
  lambda.prev <- opt.ret[['lambda1']] 
  if(opt.ret$swap) {
    lambda.prev <- opt.ret[['alpha1']] 
  }
  if(opt.ret$MAP_autofit != -1 &  iter > 1)
  {
    message(paste0("Current lambda: ", opt.ret$lambda1))
    opt.ret[['lambda1']] <- MAPfitLambda(V,opt.ret$MAP_autofit, opt.ret)
    message(paste0("Updated lambda: ", opt.ret[['lambda1']]))
    #L
    message(paste0  ("Current alpha: ", opt.ret$alpha1))
    opt.ret[['alpha1']] <- MAPfitAlpha(U,opt.ret$MAP_autofit,opt.ret)
    message(paste0("Updated alpha: ", opt.ret[['alpha1']]))
    
  }
  if(opt.ret$MAP_autofit != -1 & !is.na(opt.ret$fix.alt.setting)) #opt.ret to fix one and modulate the other, given we are doing autofit.
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
  if(length(non_empty_u) == 0){
    message('Finished');
    message('L is completely empty, alpha1 too large')
    #V = NULL;
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


CheckLowPVE <- function(X,W,V, thresh = "default") #this might be high but oh well...
{
  if(thresh == "default")
  {
    thresh = 0.005
  }else if(thresh == "avg")
  {
    thresh = 1/ncol(X)
  }
  else
  {
    thresh = 0.01
  }
  pve=PercentVarEx(as.matrix(X)*as.matrix(W), v = V)
  print(pve)
  return(which(pve <= thresh))
  
}   

#helper to just get rid of some columns quick...
DropSpecificColumns <- function(drops, mf, ms)
{
  U = as.matrix(as.data.frame(mf[, -drops]));
  V  = as.matrix(as.data.frame(ms[, -drops]));
  return(list("U"=U, "V"=V))
}

DropLowPVE <- function(X,W,V)
{
  drop.cols <- CheckLowPVE(X,W,V)
  rv <- V
  if(length(drop.cols) > 0)
  {
    rv <- V[,-drop.cols]
  }
  
  return(rv)
}

ZeroLowPVE <- function(X,W,V)
{
  drop.cols <- CheckLowPVE(X,W,V)
  message("Currently including PVE check")
  rv <- V
  rv[,drop.cols] <- 0
  return(rv)
}
#Refactor thi
AlignFactorMatrices <- function(X,W,U, V)
{

  non_empty_v = which(apply(V, 2, function(x) sum(x!=0)) > 0)
  non_empty_u = which(apply(U, 2, function(x) sum(x!=0)) > 0)
  non_empty = intersect(non_empty_u,non_empty_v);
  if(length(non_empty) < ncol(U))
  {
    message("dropping now...")
  }
  U = as.matrix(as.data.frame(U[, non_empty]));
  V  = as.matrix(as.data.frame(V[, non_empty]));
  return(list("U"=U, "V"=V))
}
# converge if: 1). Change of the values in factor matrix is small. ie. The factor matrix is stable. 2). Change in the objective function becomes small; 3). reached maximum number of iterations

ConvergenceConditionsMet <- function(iter,X,W, U,V,tracker,option)
{
  #TODO double check tracker$V tracks with the old
  #Condition 1: Change in V is small
  V_change = MatrixChange(V, tracker$V)
  if(option[['convF']] >= V_change){
    updateLog(paste0('Factor matrix converges at iteration ', iter), option);
    return(TRUE)
  }
  #option 2: objective change is small
  #Objective hasn't been updated yet in tracker. That's the issue
  
  obj_updated = compute_obj(X, W, U, V, option);
  #This isnt right
  objective_change = tracker$obj[length(tracker$obj)]- obj_updated; #newer one should be smaller than previous
  obj.change.percent <- objective_change/tracker$obj[length(tracker$obj)]
  updateLog(paste0("Objective change: ", objective_change), option)
 
  
 #If we have completed at least 1 iteration and we go up, end it.
  if(objective_change < 0 & length(tracker$obj) > 3)
  {
    message("warning: negative objective")
    print(tracker$obj)
    updateLog(paste0("Objective change going in the wrong direction (negative), ending now."), option)
    return(TRUE)
  }
  #updated change- objective change as a percent:
  if(obj.change.percent < as.numeric(option[['conv0']]) & length(tracker$obj) > 3){
    updateLog(("Objective function change threshold achieved!"), option)
    updateLog(paste0('Objective function converges at iteration ', iter), option);
    #If objective change is negative, must end....

    return(TRUE)
  }
  
  
  #Option 3- maxi number of iter.
  if(iter == option[['iter']]){
    message('Reached maximum iteration.');
    return(TRUE)
  }
  FALSE
}

MatrixChange <- function(new, old)
{
  if(ncol(new) != ncol(old))
  {
    d <- ncol(old) - ncol(new)
    nm <- matrix(rep(0, d*nrow(old)), nrow = nrow(old), ncol = d)
    new <- cbind(new,nm)
  }
  norm(new - old, 'F') / ncol(new)
}


#Main workhorse function
Update_FL <- function(X, W, option, preV = NULL, preL = NULL){
  # number of features - to avoid using T in R
  D = ncol(X)
  tStart0 = Sys.time()
  V = NULL 
  #Random initialization of F
  if(option$l_init != "")
  {
    U <- initU(X,W,option,preU=preL)
    message("in here...")
    V = fit_F(X, W, L, option)$V
      
  } else{ #initialize by F as normal.
    
    V <- initV(X,W,option,preV=preV)
  }  
  message(""); message('Start optimization ...')
  message(paste0('K = ', (option[['K']]), '; alpha1 = ', round(option[['alpha1']], digits =  4),
                 '; lambda1 = ', round(option[['lambda1']], digits = 4)))

  #Case where you use previous iteration estimates to inform next iteration..
  og_option <- option[['carry_coeffs']]
  option[['carry_coeffs']] <- FALSE
  U = FitUWrapper(X,W,V, option)
  U <- U$U
  option[['carry_coeffs']] <- og_option
  
  #Start tracking stats
  tracking.data <- UpdateTrackingParams(NULL, X,W,U,V,option)
  #If U is already empty, than we know that the matrix is too sparse, and we should just end there.
  if(CheckUEmpty(U)) {message("U is empty; ending");return(tracking.data)}
  #If we are doing a measure of per-trait variance over time...
  trait.var <- matrix(NA, option[['iter']], ncol(X))
 
  for (iii in seq(1, option[['iter']])){ #start the iterations (although technically, we have done 1 already.)
    message(paste0("Currently on iteration ", iii))
    
    ## If we are doing an autofit setting....
    if(option$MAP_autofit > -1) {
      print(iii)
      option <- UpdateSparsityMAPAutofit(iii, U,V, option)
    }
    V = fit_F(X, W, U, option, V)$V; #by iter 3 really slows down, due to the L1 requirements. Yea this won't do....
    
    #get the factor specific variance....
    if(option$traitSpecificVar)
    {
      trait.var[iii,] <- V$r.v
      V <- V$V
    }
    # if number of factors decrease because of empty factor, the change in ||F||_F = 100
    
    #Tracking change in F....
    if(CheckVEmpty(V)) {message("V is empty; ending");  return(UpdateTrackingParams(tracking.data, X,W,U,V,option))}
    colnames(V) = seq(1, ncol(V));
 
    
    ## update L
    U <- FitUWrapper(X,W,V, option,r.v = trait.var[iii,])
    if(length(U$redundant_cols) > 0)
    {
      message("Dropping some cols...")
      print(U$redundant_cols)
      r <- DropSpecificColumns(U$redundant_cols, U$U, V)
      U <- r$U; V <- r$V
    }
      else{
        U <- U$U
      }
  
    # if L is empty, stop
    if(CheckUEmpty(U)) {message("U is empty; ending"); return(UpdateTrackingParams(tracking.data, X,W,U,V,option))} #Update and end.
    colnames(U) = seq(1, ncol(U)) ;
    
    # Drop low PVE and align two matrices	
    #V <- DropLowPVE(X,W,V)
    updated.mats <- AlignFactorMatrices(X,W,U, V); U <- updated.mats$U; V <- updated.mats$V


    if(ConvergenceConditionsMet(iii,X,W,U, V, tracking.data, option))
    {
      message("Convergence criteria met...")
      
      #Check out low pve one..
      drops <- CheckLowPVE(X,W,V) #redundant with other code
      print(drops)
      if(length(drops) > 0)
      {
        message("Pruning out low PVE columns")
        n <- DropSpecificColumns(drops, U,V)
        U <- n$U; V <- n$V
      }
      cat('\n')
      updateLog(paste0('Total time used for optimization: ',  round(difftime(Sys.time(), tStart0, units = "mins"), digits = 3), ' min'), option);
      cat('\n')
     return(UpdateTrackingParams(tracking.data, X,W,U,V,
                                            option))
    }else{
      tracking.data <- UpdateTrackingParams(tracking.data, X,W,U,V,
                                            option)
      
      if(option$V > 0){
        EndIterStatement(iii, tracking.data, option)
        }
    }
  }
}

EndIterStatement <- function(iter, td, option)
{
  cat('\n')
  #Update message- need to change this....
  updateLog(paste0('Iter', iter, ':'), option)
  updateLog(paste0('Objective change = ', td$obj_change[length(td$obj_change)]), option)
  updateLog(paste0('Frobenius norm of (updated factor matrix - previous factor matrix) / number of factors  = ', td$V_change[length(td$V_change)]), option);
  updateLog(paste0('Loading Sparsity = ', td$U_sparsities[length(td$U_sparsities)], '; Factor sparsity = ', td$V_sparsities[length(td$V_sparsities)], '; ', td$K, ' factors remain'), option); 
  cat('\n')
}


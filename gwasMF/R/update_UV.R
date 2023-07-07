################################################################################################################################
#Modified 4/26 by Ashton Omdahl
## Copied over from another function since this wasn't working originally.
################################################################################################################################

#helpful code for debugging
if(FALSE)
{
  #source("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/quickLoadData.R")
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/quickLoadData.R")
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_V.R")
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
1/(NROW(x)-1) * crossprod(scale(x, TRUE, TRUE))
}

#This function uses OLS to get a good estimate of what the maximum sparsity space is across all parameters.
#@param burn.in- how many iterations to go
#@return V_burn: the burn-in estimate of V
#@return U_burn: the burn-in estimate of U
#@return max_sparsity_params: a list indexed by iteration containing the list of max alphas and lambdas (i.e. list[[i]]$alpha))
#DefineSparsitySpaceInit(X, W, option, burn.in = burn.in.iter)
DefineSparsitySpaceInit <- function(X, W,W_c, W_ld, option, burn.in = 5,...)
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
  V.dat <- initV(Xint,Wint, new.options, ...)
  V<-V.dat$V

  if(burn.in < 1)
  {
    if(option$u_init != "")
    {
      message("initializing with U")
      U <- initU(Xint, Wint, new.options)
      V.dat <- DefineSparsitySpace(Xint, Wint, W_c, U, "V", new.options, fit = "None")
      new.options$K = ncol(U)
      param.space[[1]] <- list("alpha" = NULL, "lambda"= V.dat)
      return(list("V_burn" =NULL , "U_burn"=U, "max_sparsity_params"=param.space, "new.k" =new.options$K))
    }
    else
    {
      #2-24 changes
      #if(option$fixed_ubiq)
      #{
      #  U.dat <- DefineSparsitySpace(Xint, Wint, V[,-1], "U", new.options, fit = "None") #this is associated with alpha
      #}else
      #{
        U.dat <- DefineSparsitySpace(Xint, Wint,W_c, V, "U", new.options, fit = "None")
      #}
      new.options$K = ncol(V)
      param.space[[1]] <- list("alpha" = U.dat, "lambda"= NULL)
      return(list("V_burn" = V, "U_burn"=NULL, "max_sparsity_params"=param.space, "new.k" =new.options$K))
    }


    #V.dat <- DefineSparsitySpace(Xint, Wint,U.dat$U, "V", new.options, fit = "none")
    #param.space[[1]] <- list("alpha" = U.dat$sparsity_space, "lambda"= V.dat$sparsity_space)
    #V.dat$sparsity_space)
    #return(list("V_burn" = V, "U_burn"=U.dat$U, "max_sparsity_params"=param.space, "new.k" =new.options$K))

  }
  for(i in 1:burn.in)
  {
    print(i)
    #U.dat <- FitUWrapper(Xint,Wint,V, new.options)
    U.dat <- DefineSparsitySpace(Xint, Wint, W_c, V, "U", new.options, fit = "OLS")
    new.options$K = ncol(U.dat$U) #Update if it has changed
    #If we have ones with NA, they need to get dropped
    #if we have columns with NAs here, we want them gone now so it doesn't jank up downstream stuff.
    #V.dat <- FitVWrapper(Xint, Wint, U.dat$U, new.options, V);
    V.dat <- DefineSparsitySpace(Xint, Wint,W_c,U.dat$U, "V", new.options, fit = "OLS")
    if(i > 3)
    {
      #this might be high but oh well...
      #HERE
      #drops <- CheckLowPVE(X,W,V.dat$V) #redundant with other code
      drops <- c()
      if(length(drops) > 0)
      {
        n <- DropSpecificColumns(drops, U.dat$U, V.dat$V)
        U.dat$U <- n$U
        new.options$K <- ncol(U.dat$U)
        V.dat <- DefineSparsitySpace(Xint, Wint,W_c,U.dat$U, "V", new.options, fit = "OLS")
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
#TODO: makea. test here to check the dimensions of W_c and X and fixed. Needs to line upo.
#' Title
#'
#' @param X
#' @param W
#' @param W_cov
#' @param fixed
#' @param learning
#' @param option
#' @param fit
#'
#' @return
#' @export

DefineSparsitySpace <- function(X,W,W_cov,fixed,learning, option, fit = "None")
{
  new.options <- option
  new.options$regression_method <- fit
  new.options$actively_calibrating_sparsity <- TRUE

  if(learning == "V")
  {
    #Not implemented W_cov adjustment for V at this time...
    free.dat <- FitVWrapper(X, W, W_cov, fixed, new.options);
    #HERE
  }else #learning U
  {
    #Check: ensure
    stopifnot(!is.null(W_cov))
    stopifnot(nrow(W_cov) == ncol(X))
    red.cols <- c(1)
    while(!is.null(red.cols))
    {

      free.dat <- FitUWrapper(X,W, W_cov, fixed, new.options)
      red.cols <- free.dat$redundant_cols
      if(!is.null(red.cols)){
        message("capturing redundant columns...")
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
initV <- function(X,W,option, preV = NULL, rg_ref = NULL)
{
  D = ncol(X)
  cor_struct <- cor2(X)
  #svd.r <- svd(cor_struct, nu = option$K)
  svd.corr <- svd(cor_struct)
  svd.dat <- svd(X)
  message("Wastefully calculating all SVs now, change this later.")
  setK = selectInitK(option,X,W,svs = svd.r$d)
  if(!option[['preinitialize']])
  {
    message("TODO: initializing only positive for now.")
    V = matrix(runif(D*(setK - 1), min = -1, max = 1), nrow = D);

    if(option[['f_init']] == 'ones_plain')
    {
      message("Initializing ubiquitous factor with all ones...")
      V = cbind(rep(1, D), V);
    }	else if(option[['f_init']] == 'ones_eigenvect') {
      message("1st column based on direction of svd of cor")
      ones <- sign(svd.corr$u[,1])
      #Reintroduced later- no problems here, and gives a nice speedup I think.
      if(option$svd_init)
      {
        message("initializing with SVD")
        V <- svd.dat$v[,1:(setK-1)]
        #V <- svd$u[,2:(option$K)]
      } #otherwise its random.

    } else if(option[['f_init']] == 'plieotropy')
    {
      message("1st column based svd of cor(|Z|), since plieotropy has no direction.")
      ones <- svd.r$u
    }else if(grepl(pattern = "rg_ref",x = option$f_init) & !is.null(rg_ref)) #Use an RG reference file
    {
      message("using the passed in rg...")
      evals <- base::eigen(rg_ref)
      ones = evals$vectors[,1]
      V <- evals$vectors[,2:(setK)]
      if(grepl(pattern = 'ones_eigenvect',x = option$f_init))
      {
        ones <- sign(svd.r$u[,1])
        V <- evals$vectors[,1:(setK-1)]
      }

    }
    else {
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
  s=1
  return(list("V"=V, "s"=s))
}

#Function to inialize U if specified
initU <- function(X,W,option, prevU = NULL)
{
  message("Initializing L rather than F.")
  nsnps = nrow(X)
  U   = matrix(runif(nsnps*(option[['K']] - 1)), nrow = nsnps);
  if(option$u_init == 'ones_plain')
  {
    message("Initializing ubiquitous factor with all ones...")
    ones = rep(1, nsnps)

  }	else if(option$u_init == 'ones_eigenvect') {
    message("1st column based on direction of svd of cor")
    cor_struct <- cor2(t(X))
    svd <- irlba::irlba(cor_struct, 1) #fortunately its symmetric, so  U and V are the same here!
    ones <- sign(svd$u)
  } else if(option$u_init == 'plieotropy') {
    message("1st column based svd of cor(|Z|), since plieotropy has no direction.")
    cor_struct <- cor2(abs(t(X)))
    svd <- irlba(cor_struct, 1) #fortunately its symmetric, so  U and V are the same here!
    ones <- svd$u
  } else {
    ones = matrix(runif(nsnps*(option[['K']])), nrow = nsnps)[,1]
  }
  if(option$svd_init)
  {
    intermediate <- svd(X, nu = (option[['K']]))
    U <- intermediate$u[,2:(option[['K']])]
    ones <- intermediate$u[,1]
    message("Initializing U with a fit from SVD")
  }
  cbind(ones, U);

}

#helper code to clean things up, just keep track of those relevant metrics,
UpdateTrackingParams <- function(sto.obj, X,W,W_c,U,V,option, sparsity.thresh = 1e-5, loglik = NULL, scalar=1)
{
  if(is.null(sto.obj))
  {
    sto.obj <-  list("V" = V, "U" = U, "initK" = option$K, "K" = ncol(U), "obj" = c(NA), "obj_change" = c(),
                     "V_sparsities" = c(), "U_sparsities" = c(), "autofit_lambda" = c(), "autofit_alpha"=c(), "mse" = c(),
                     "V_change" = c(), "U_change" = c(), "decomp_obj" ="", "model.loglik" = c(), "Vs"=list(), "Us"=list())
  }
    # collect sparsity in L and F

  #sto.obj$U_sparsities = c(sto.obj$U_sparsities, sum(abs(U) < sparsity.thresh) / (ncol(U) * nrow(U)));
  #sto.obj$V_sparsities = c(sto.obj$V_sparsities, sum(abs(V) < sparsity.thresh) / (ncol(V) * nrow(V)));
  sto.obj$U_sparsities = c(sto.obj$U_sparsities, matrixSparsity(U, option$K)); #The sparsity count should be with respect to the initial value..
  #TODO: fix the relative sparsity count thuing...
  sto.obj$V_sparsities = c(sto.obj$V_sparsities,matrixSparsity(V, option$K))
  sto.obj$K = ncol(U);
  sto.obj$mse <- c(sto.obj$mse, norm(t(W_c %*% t(W*X))-t(W_c %*% t(W*(U%*%t(V)))), type = "F")/(nrow(X) * ncol(X)))
    # change in the objective -- should always be a decrease
  #Objective and fit and lll
  if(!is.null(loglik))
  {
    sto.obj$model.loglik <- c(sto.obj$model.loglik,loglik)
  }
    obj_updated = compute_obj(X, W,W_c, U, V, option, scalar=scalar)
    sto.obj$decomp_obj = compute_obj(X, W,W_c, U, V, option,decomp = TRUE, scalar=scalar)
  if(is.na(sto.obj$obj[1]))
  {
    sto.obj$obj[1] <- obj_updated
  }
  obj_change = sto.obj$obj[length(sto.obj$obj)] - obj_updated;
  sto.obj$obj = c(sto.obj$obj, obj_updated);
  sto.obj$obj_change = c(sto.obj$obj_change, obj_change);
    #Update the fit:

  sto.obj$V_change = c(sto.obj$V_change, MatrixChange(V, sto.obj$V))
  sto.obj$U_change = c(sto.obj$U_change, MatrixChange(U, sto.obj$U))
  #Update the sparsity paramters, regardless of if they change or not.
      sto.obj$autofit_lambda <- c(sto.obj$autofit_lambda,option$lambda1)
      sto.obj$autofit_alpha <- c(sto.obj$autofit_alpha,option$alpha1)

      #LAST STEP- update the new U and V
      sto.obj$U <- as.matrix(U)
      sto.obj$Us[[length(sto.obj$Us) + 1]] <-  as.matrix(U)
      sto.obj$V <- as.matrix(V)
      sto.obj$Vs[[length(sto.obj$Vs) + 1]] <-  as.matrix(V)
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

DropLowPVE <- function(X,W,V,...)
{
  drop.cols <- CheckLowPVE(X,W,V,...)
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
#If the objective goes negative, you want the old one
ConvergenceConditionsMet <- function(iter,X,W,W_c, U,V,tracker,option, initV = FALSE, loglik = NULL, scalar = 1)
{
  #TODO double check tracker$V tracks with the old
  #1 conveged
  #2 exceeded
  #3 not_converged
  conv.opts <- c("converged", "exceeded", "not_converged")
  #Condition 1: Change in V is small
  V_change = MatrixChange(V, tracker$V)
  if(option[['convF']] >= V_change){
    updateLog(paste0('Factor matrix converges at iteration ', iter), option);
    return(conv.opts[[1]])
  }
  #option 2: objective change is small
  #Objective hasn't been updated yet in tracker. That's the issue

  obj_updated = compute_obj(X, W,W_c, U, V, option, scalar=scalar)
  #This isnt right
  objective_change = tracker$obj[length(tracker$obj)]- obj_updated; #newer one should be smaller than previous
  obj.change.percent <- objective_change/abs(tracker$obj[length(tracker$obj)])
  #Special case- we initialized with a V, and our next iteration is worse
  if(objective_change < 0 & initV & length(tracker$obj) == 1)
  {
    message("Previous iterations converged already. Ending now.")
    print(c(tracker$obj, obj_updated))
    updateLog(paste0("Objective change going in the wrong direction (negative), ending now."), option)
    return(conv.opts[[2]])
  }
 #If we have completed at least 1 iteration and we go up, end it.
  #First iteration has objective change of 0.
  if(objective_change < 0 & option[['conv0']] > 0) #& length(tracker$obj) > 2)
  {
    message("warning: negative objective")
    print(tracker$obj)
    updateLog(paste0("Objective change going in the wrong direction (negative), ending now."), option)
    return(conv.opts[[2]])
  }
  #updated change- objective change as a percent:
  #TODO: add condition to see if it has a prior on it or not. This will shorten our number of iterations.
  if((obj.change.percent <= as.numeric(option[['conv0']])) & (length(tracker$obj) > 2)){
    updateLog(("Objective function change threshold achieved!"), option)
    updateLog(paste0('Objective function converges at iteration ', iter), option);
    #If objective change is negative, must end....
    return(conv.opts[[1]])
  }


  #Option 3- maxi number of iter.
  if(iter == option[['iter']]){
    message('Reached maximum iteration.');
    return(conv.opts[[1]])
  }
  conv.opts[[3]]
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
#1.26 changes- drop the burn in.
Update_FL <- function(X, W, W_c, option, preV = NULL, preU = NULL, burn.in = FALSE){
  # number of features - to avoid using T in R
  #Tracking data for debugging things...
  ll.tracker <- c()
  mse.tracker <- c()
  penalty.tracker <- c()
  iter.tracker <- c()
  D = ncol(X)
  tStart0 = Sys.time()
  es.objective <- c()
  em.objective <- list()
  iter.objective <- list()
  s.weight.tracker <- list()
  easy.objective <- c()
  dev.score<- c()
  V = NULL
  #Random initialization of F
  if(!is.null(preV))
  {
    V = preV
  }
  else if(option$u_init != "")
  {
    U <- initU(X,W,option,preU=preU)
    message("initializing with U")
    V.dat =  fit_V(X, W, L, option)
    V <- V.dat$V/V.dat$s

  } else if(burn.in)
  {
    burn.in.sparsity <- DefineSparsitySpaceInit(X, W,W_c,NULL, option, burn.in = 4)
    V <- burn.in.sparsity$V_burn
  }
  else{ #initialize by F as normal.

    V.dat <- initV(X,W,option,preV=preV) #returns scalar and V at unit norm length, #CONFIRM
    V <- V$V
  }
  message(""); message('Start optimization ...')
  message(paste0('K = ', (option[['K']]), '; alpha1 = ', round(option[['alpha1']], digits =  4),
                 '; lambda1 = ', round(option[['lambda1']], digits = 4)))

  #Case where you use previous iteration estimates to inform next iteration..
  og_option <- option[['carry_coeffs']]
  option[['carry_coeffs']] <- FALSE

  U.dat = FitUWrapper(X,W,W_c,V, option) #returns it scaled, with S  #CONFIRM
  U <- U.dat$U
  option[['carry_coeffs']] <- og_option
  #if(option$debug)
  #{em.objective[[0]] <- compute_obj(X,W,W_c, U,V,option, globalLL=TRUE,decomp = TRUE)}
  #Start tracking stats
  tracking.data <- UpdateTrackingParams(NULL, X,W,W_c,U,V,option,scalar=U.dat$s) #multiplies everything by the scalar: #CONFIRM

  #If U is already empty, than we know that the matrix is too sparse, and we should just end there.
  if(CheckUEmpty(U)) {message("U is empty; ending");return(tracking.data)}

  #If we are doing a measure of per-trait variance over time...
  trait.var <- matrix(NA, option[['iter']], ncol(X))
  i.c <- 1
  for (iii in seq(1, option[['iter']])){ #start the iterations (although technically, we have done 1 already.)
    message(paste0("Currently on iteration ", iii))
    iteration.ll.total <- 0
    ## If we are doing an autofit setting....
    if(option$MAP_autofit > -1) {
      print(iii)
      option <- UpdateSparsityMAPAutofit(iii, U,V, option)
    }

    V.new = FitVWrapper(X, W,W_c, U, option)#, formerV = V);  #returns V scaled, with S  #CONFIRM
    iteration.ll.total <- V.new$total.log.lik
    ll.tracker <- c(ll.tracker, V.new$total.log.lik)
    penalty.tracker <- c(penalty.tracker, V.new$penalty)
    V.prev <- V
    V = V.new$V
    #The scaling term applied to U
    mse.tracker <- c(mse.tracker, norm(W*X-W*(U%*%t(V))/V.new$s, type = "F")/(nrow(X) * ncol(X)))
    iter.tracker <- c(iter.tracker, "V")
    easy.objective <- c(easy.objective, getQuickObj(V.new$SSE,nrow(X) * ncol(X), U,V,option))
    dev.score <- c(dev.score, V.new$SSE)
    #More crap to follow the objective
    #if(option$debug)
    #{
      es.objective <- c(es.objective, GetStepWiseObjective(X,W,W_c,V.prev,V*V.new$s, U,"U",option));
      em.objective[[(i.c)]] <- compute_obj(X,W,W_c, U,V,option, globalLL=TRUE, decomp = TRUE, scalar=V.new$s)
      s.weight.tracker[[i.c]] <- V.new$s
      i.c <- i.c + 1
    #}
    #Just get the matrix out.
    #get the factor specific variance....
    if(option$traitSpecificVar)
    {
      trait.var[iii,] <- V.new$r.v
    }
    # if number of factors decrease because of empty factor, the change in ||F||_F = 100

    #Tracking change in F....
    if(CheckVEmpty(V)) {message("V is empty; ending");  return(UpdateTrackingParams(tracking.data, X,W,W_c,U,V,option, loglik = iteration.ll.total,scalar=V.dat$s))}
    colnames(V) = seq(1, ncol(V));


    ## update L
    U.new <- FitUWrapper(X,W,W_c,V, option,r.v = trait.var[iii,])#, prevU = U) #returns U scaled, with S  #CONFIRM
    U.prev <- U
    U = U.new$U #confirm- its unit norm?
    ll.tracker <- c(ll.tracker, U.new$total.log.lik)
    penalty.tracker <- c(penalty.tracker, U.new$penalty)
    mse.tracker <- c(mse.tracker, norm(W*X-W*(U %*%t(V))/U.new$s, type = "F")/(nrow(X) * ncol(X)))
    easy.objective <- c(easy.objective, getQuickObj(U.new$SSE,nrow(X) * ncol(X), U,V,option))
    dev.score <- c(dev.score, U.new$SSE)
    iter.tracker <- c(iter.tracker, "U")
    #if(option$debug)
    #{
      #ES stands for each step.
      es.objective <- c(es.objective, GetStepWiseObjective(X,W,W_c,U.prev,U*U.new$s, V,"V",option))
      em.objective[[(i.c)]] <- compute_obj(X,W,W_c, U,V,option, globalLL=TRUE, decomp = TRUE, scalar= U.new$s)
      iter.objective[[(i.c/2)]] <- compute_obj(X,W,W_c, U,V,option, globalLL=TRUE, decomp = TRUE, scalar=U.new$s)
      s.weight.tracker[[i.c]] <- U.new$s
      i.c <- i.c + 1
    #}

    iteration.ll.total <- iteration.ll.total + U.new$total.log.lik
    if(length(U.new$redundant_cols) > 0)
    {
      message("Dropping some cols...")
      print(U.new$redundant_cols)
      r <- DropSpecificColumns(U.new$redundant_cols, U, V)
      U <- r$U; V <- r$V
    }
    # if L is empty, stop
    #rescale V, since we last learned U with a rescaled V
    if(CheckUEmpty(U)) {message("U is empty; ending"); return(UpdateTrackingParams(tracking.data, X,W,W_c,U,V,option,loglik = iteration.ll.total,scalar=U.dat$s))} #Multiplies by the extra factor.
    colnames(U) = seq(1, ncol(U)) ;

    # Drop low PVE and align two matrices
    #V <- DropLowPVE(X,W,V)
    updated.mats <- AlignFactorMatrices(X,W,U, V); U <- updated.mats$U; V <- updated.mats$V

    convergence.status <- ConvergenceConditionsMet(iii,X,W,W_c,U, V,
                                                   tracking.data, option, initV = is.null(preV),
                                                   loglik = iteration.ll.total, scalar = U.dat$s)
    if(convergence.status %in% c("converged", "exceeded"))
    {
      message("Convergence criteria met...")
      cat('\n')
      updateLog(paste0('Total time used for optimization: ',  round(difftime(Sys.time(), tStart0, units = "mins"), digits = 3), ' min'), option);
      cat('\n')
      if(convergence.status == "exceeded")
      {
        #return the previous iteration, not the next
        ret <- UpdateTrackingParams(tracking.data, X,W,W_c,tracking.data$U,tracking.data$V,
                                    option, loglik = iteration.ll.total, scalar = U.dat$s) #CONFIRM
      }else
      {
        ret <- UpdateTrackingParams(tracking.data, X,W,W_c,U,V,
                                    option, loglik = iteration.ll.total,scalar=U.dat$s) #CONFIRM
      }
      #if(option$debug){
        ret$each.step.obj <- es.objective;
        ret$each.matrix.obj <- em.objective;
        ret$full.iter.obj <- iter.objective
        ret$penalized_ll <- ll.tracker
        ret$penalties <- penalty.tracker
        ret$iter_mat <- iter.tracker
        ret$global.mse <- mse.tracker
        ret$final.scaling <- U.new$s
        ret$dev.obj <- easy.objective
        ret$dev.only <- dev.score
        ret$s_V <- U.new$s #this was the scalar used to adjust V, in the last step
        ret$s_U <- V.new$s #This was the scalar used to adjust U
       # }
     return(ret)
    }else{
      tracking.data <- UpdateTrackingParams(tracking.data, X,W,W_c,U,V,
                                            option, loglik = iteration.ll.total,scalar=U.dat$s)

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
  updateLog(paste0('Percent objective change = ', abs(td$obj[length(td$obj)-1]-td$obj[length(td$obj)])/td$obj[length(td$obj)]), option)
  updateLog(paste0('Frobenius norm of (updated factor matrix - previous factor matrix) / number of factors  = ', td$V_change[length(td$V_change)]), option);
  updateLog(paste0('U Sparsity = ', round(td$U_sparsities[length(td$U_sparsities)], digits = 3),
                   '; V sparsity = ',  round(td$V_sparsities[length(td$V_sparsities)],digits = 3), '; ', td$K, ' factors remain'), option);
  cat('\n')
}


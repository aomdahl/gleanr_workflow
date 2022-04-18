  ################################################################################################################################
  ## Update the factor matrix L during the alternative optimization
  ##
  ## Input:
  ##              - X: the matrix to be decomposed (N x T, N is the number of data points, T is the number of features)
  ##              - F: learned factor matrix  (T x K, T is the number of features, K is the number of factors)
  ##              - W: the weight matrix, same size as in X
  ##              - option: a list with parameters including lambda1, the l1 penalty parameter for the factor matrix (F)
  ##
  ## Return:
  ##              - A loading matrix with mixed signs that minimize_F ||(X - LF') .* W||_F^2 + lambda1*|L|_1
  ##
  ## Example of usage:
  ##
  ## source('../simulation/Generate_input.R');
  ## data = generate_input(tau = tau);
  ## F = data[[1]];
  ## X = data[[3]];
  ## W = data[[4]];
  ## option = list();
  ## option[['lambda1']] = 0.1;
  ## L_predict = fit_L(X, W, F, option);
  ##
  ################################################################################################################################
  
  
  suppressWarnings(library(penalized))
  library(foreach)
  library(parallel)
  
  #Function for running the fitting step. Options include
  #reweighted: this uses the previous iteration to initialize the current one. Thought it might provide a speed boost, but the difference seemed to be minimal.
  #glmnet: This uses the glmnet function; just exploring options
  #fastReg: This one was supposed to be faster. I don't know that it was at all; but not sure about serious results on it
  #ridge_L: this employs L2 as opposed to L1 regression  on L
  #fixed_ubiq: this allows for a ubiquitous factor by removing the L1 constraint from the first column.
  #@param r.v- the trait-specific variance vector.
  #@return: l, the current learned row of l (k x 1) of weights for a given SNP
  one_fit <- function(x, w, FactorM, option, formerL, r.v){
    if(option$traitSpecificVar && !is.null(r.v))
    {
      xp = x / sqrt(r.v);
      #print("In here...")
      #xp = x #I don't think we weight the samples by the variance, do we?
      FactorMp = diag(1/sqrt(r.v)) %*% FactorM;
    }else{
      xp = unlist(w * x);
      FactorMp = diag(w) %*% FactorM; #scaling this way is much smaller
    }
     #elementwise multiplication
    #if(option$parallel)
    #{
    #  xp = t(w * x)
    #}

    # Fit: xp' = FactorMp %*% l with l1 penalty on l -- |alpha1 * l|
    #Removed the transpose on X, getting the wrong dimensions
    dat_i = as.data.frame(cbind((xp), FactorMp));
    colnames(dat_i) = c('X', paste0('F', seq(1, ncol(FactorMp))));

    if(option[["reweighted"]])
    {
      fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=1e-10,
                      positive = FALSE, standardize = FALSE, trace = FALSE, startbeta = formerL);
      l = coef(fit, 'all');
    } 
    else if(option[["glmnet"]]) 
      {
      
      fit <- glmnet(x = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], y = xp, alpha = 1, lambda = c(option[['alpha1']]), intercept = FALSE)
      l = coef(fit, 'all')[,1][-1];
      
    }	
    else if(option[["fastReg"]]){
      td <- t(dat_i[,paste0('F', seq(1,ncol(FactorMp)))])
      XTX <- td %*% dat_i[,paste0('F', seq(1,ncol(FactorMp)))]
      XTY <- td %*% xp
      fit <- elasticnet(XTX, XTY, lam2 = -1, lam1 = option[['alpha1']])
      #did we really try this?
      
    }	
    else if(option[["ridge_L"]]){
      
      print("doing the ridge L")
      fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 = 1e-15, lambda2=option[['alpha1']],
                      positive = FALSE, standardize = FALSE, trace = FALSE);
      l = coef(fit, 'all');
      
    } 
    else if(option$calibrate_sparsity)
    {
      fit = recommendRange(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(FactorMp)))], data=dat_i,
                           unpenalized = ~0, lambda1 =lambdas, lambda2=1e-10,
                           positive = FALSE)
      return(fit)
    }
    else if(option[["fixed_ubiq"]])
    {
      lambdas <- c(0, rep(option[['alpha1']], (ncol(FactorMp) - 1))) #no lasso on that first column
        #browser()
      #Okay, it seems like this is working now to estimate the max. Track over all the rows, then get the max
      fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 =lambdas, lambda2=1e-10,
                      positive = FALSE, standardize = FALSE, trace = FALSE)
      l = coef(fit, 'all')
    }		else {
      fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=1e-10,
                      positive = FALSE, standardize = FALSE, trace = FALSE);
      l = coef(fit, 'all');
    }
    
    return(l)
  }
  
  
  fit_L<- function(X, W, FactorM, option, formerL, r.v = NULL){
  	L = NULL
  	tS = Sys.time()
  	for(row in seq(1,nrow(X))){
  		x = X[row, ];
  		w = W[row, ];
  		if(option[['reweighted']])
  		{
  		  l = one_fit(x, w, as.matrix(FactorM), option, formerL[row,], r.v);
  		}		else {
  		  l = one_fit(x, w, as.matrix(FactorM), option, NULL, r.v);
  		}
  		L = rbind(L, l);
  		#print(row)
  	}
  	tE = Sys.time();
  	print(paste0('Updating Loading matrix takes ', round((tE - tS)/60, 2), 'min'));
  
  	return(L)
  }
  

#Something is jacked up with this.
  fit_L_parallel <- function(X, W, FactorM, option, formerL){
    L = NULL
    tS = Sys.time()
    cl <- parallel::makeCluster(option[["ncores"]])
    print("c-1")
    doParallel::registerDoParallel(cl)
    L <- foreach(row =seq(1,nrow(X)), .combine = 'rbind', .packages = 'penalized') %dopar% {
      x = X[row, ];
      w = W[row, ];
      xp = matrix(w * x, nrow = nrow(FactorM), ncol = 1); #elementwise multiplication
      FactorMp = diag(w) %*% FactorM;  #what are we doing here?
      dat_i = data.frame(cbind(xp, FactorMp));
      colnames(dat_i) = c('X', paste0('F', seq(1, ncol(FactorMp))));
      fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=1e-10,
                      positive = FALSE, standardize = FALSE, trace = FALSE);
      l = coef(fit, 'all');
      #should be able to do with a one_fit command here, once we figure out why this isn't working....
      l
    }
    
    tE = Sys.time();
    message(paste0('Updating Loading matrix takes ', round((tE - tS)/60, 2), 'min'));
    return(L)
  }

  

  ################################################################################################################################
  ## Update the factor matrix L during the alternative optimization
  ## Based originally on Yuan He's code, with some speed-up improvements and tweaks made by Ashton 2019-2022
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
  library(doParallel)
  library(svMisc)
  #Function for running the fitting step. Options include
  #carry_coeffs: this uses the previous iteration to initialize the current one. Thought it might provide a speed boost, but the difference seemed to be minimal.
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
      FactorMp = diag(1/sqrt(r.v)) %*% FactorM;
    }else{
      xp = unlist(w * x);
      FactorMp = diag(w) %*% FactorM; #scaling this way is much smaller
    }
    #if we are doing the burn in....
    max.sparsity <- NULL
    # Fit: xp' = FactorMp %*% l with l1 penalty on l -- |alpha1 * l|
    #Removed the transpose on X, getting the wrong dimensions
    dat_i = as.data.frame(cbind((xp), FactorMp));
    colnames(dat_i) = c('X', paste0('F', seq(1, ncol(FactorMp))));

    if(option[["carry_coeffs"]])
    {
      fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=1e-10,
                      positive = FALSE, standardize = FALSE, trace = FALSE, startbeta = formerL);
      l = coef(fit, 'all');
    }else if (option[["regression_method"]] == "OLS") {
      fit <- lm(xp~FactorMp + 0) 
      l = coef(fit) #check this, but I think its correct
      if(option$actively_calibrating_sparsity) { max.sparsity <- rowiseMaxSparsity(xp, FactorMp)}
      
    }
    else if(option[["regression_method"]] == "glmnet" & option[["fixed_ubiq"]]) 
      {
      penalties <- c(0, rep(1, (ncol(FactorMp) - 1)))
      #fit <- glmnet(x = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], y = xp, alpha = 1, lambda = option[['alpha1']], 
                    #intercept = FALSE, penalty.factor = penalties)
      #l = coef(fit, 'all')[,1][-1];
      l = glmnetLASSO(dat_i, xp, ncol(FactorMp), option[['alpha1']], penalties)
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
    else if( option[["regression_method"]] == "penalized" & option[["fixed_ubiq"]])
    {
      #lambdas <- c(0, rep(option[['alpha1']], (ncol(FactorMp) - 1))) #no lasso on that first column
      #Above was not correct- it was supposed to be sparse,, not sure why this is what it is
      lambdas <-  rep(option[['alpha1']], ncol(FactorMp))
      fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 =lambdas, lambda2=1e-10,
                      positive = FALSE, standardize = FALSE, trace = FALSE, epsilon = option$epsilon, maxiter = 10) #seet a limit on maxiter.
      #Note- some rudimentary testing done, setting maxiter doesn't necessarily help. 
      l = coef(fit, 'all')
    }		else {
      fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=1e-10,
                      positive = FALSE, standardize = FALSE, trace = FALSE);
      l = coef(fit, 'all');
    }
    return(list("l" = l, "max_sparse"=max.sparsity))
   
  }
  
  
  fit_L<- function(X, W, FactorM, option, formerL, r.v = NULL){
  	L = NULL
  	sparsity.est <- c()
  	tS = Sys.time()
    #updateLog("Updating L...", option)
  	for(row in seq(1,nrow(X))){
      if(option$V > 0)
      {
        svMisc::progress(row, nrow(X), progress.bar = TRUE)
      }
  		x = X[row, ];
  		w = W[row, ];
  		if(option[['carry_coeffs']])
  		{
  		  l = one_fit(x, w, as.matrix(FactorM), option, formerL[row,], r.v);
  		}		else {
  		  l = one_fit(x, w, as.matrix(FactorM), option, NULL, r.v);
  		}
  		#Gettings some very bizarre errors about names and such, some weird hacks to work around it.
  		if(length(l$l) == 1)
  		{
  		  l <- as.matrix(l$l)
  		  suppressMessages(names(l) <- "F1")
  		}
      if(is.null(L))
      {
          L = l$l
      } else if(is.null(L) & is.null(l)){
          message("Matrix is empty...")
          return(NULL)
      }
      else{
          L = suppressMessages(bind_rows(L, l$l))
      }
  		if(option$actively_calibrating_sparsity) #we are in the process of burning in with OLS and calibrating sparsity
  		{
  		  sparsity.est <- c(sparsity.est, l$max_sparse)
  		}
      
  	}
  	updateLog(paste0('Updating Loading matrix takes ', round(difftime(Sys.time(), tS, units = "mins"), digits = 3), ' min'), option);
  
  	return(list("U" = as.matrix(L), "sparsity_space"=sparsity.est))
  }
  

#Should be working fine... need to do some testing though!
  fitUParallel <- function(X, W, FactorM, option, formerL){

    L = NULL
  	message("inside the parallele list")
    tS = Sys.time()
    #cl <- parallel::makeCluster(option[["ncores"]])#, outfile = "TEST.txt")
    iterations <- nrow(X)
    #doParallel::registerDoParallel(cl)
    doParallel:registerDoParallel(parallel::makeCluster(option[["ncores"]]))
    #6/13/2022
    #We split the tasks across cores, so multiple regression steps per core; you aren't switching nearly as much.
    n.per.group <- ceiling(nrow(X)/option$nsplits)
    oo <- n.per.group - 1
    print(option$nsplits)
    split_lines <- lapply(1:(option$nsplits-1), function(x) (x*n.per.group - oo) :(x*n.per.group))
    split_lines[[option$nsplits]] <- (split_lines[[(option$nsplits-1)]][n.per.group]+1):nrow(X)
    L <- foreach(rows = split_lines, .combine = 'bind_rows', .packages = c('penalized', 'dplyr')) %dopar% {
      sub_l <- NULL
      for(row in rows)
      {
        x = X[row, ];
        w = W[row, ];
        
        xp = matrix(w * x, nrow = nrow(FactorM), ncol = 1); #elementwise multiplication
        FactorMp = diag(w) %*% FactorM;  #what are we doing here?
        dat_i = data.frame(cbind(xp, FactorMp));
        colnames(dat_i) = c('X', paste0('F', seq(1, ncol(FactorMp))));
        fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                        unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=1e-10,
                        positive = FALSE, standardize = FALSE, trace = FALSE,epsilon = option$epsilon);
        sub_l <- bind_rows(sub_l, coef(fit, 'all'))
      }
      sub_l
    }
    doParallel::stopImplicitCluster()
   #parallel::stopCluster(cl) 
    updateLog(paste0('Updating Loading matrix takes ',  as.character(round(difftime(Sys.time(), tS, units = "mins"), digits = 3), ' min')), option);
    return(as.matrix(L))
  }


#      if(row %% 100 == 0)
#      {
#        sink("log.txt", append=TRUE)
#        cat(paste("Starting iteration",row,"\n"))
#        sink()
#      }

#doing some major renaming
  fit_L_parallel <- fitUParallel
  FitUWrapper <- function(X,W,FactorM, option,...)
  {
    if(option[['ncores']] > 1) #This is not working at all. Can't tell you why. But its not. Need to spend some time debugging at some point.
    {
      log_print("Fitting L in parallel")	
      U = fitUParallel(X, W, FactorM, option, formerL = preL); #preL is by default Null, unless yo specify!
      
    }else
    {
      U = fit_L(X, W, FactorM, option, ...);
    }
    return(U)
  }


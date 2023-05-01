  ################################################################################################################################
  ## Update the factor matrix U during the alternative optimization
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
  ## L_predict = fit_U(X, W, F, option);
  ##
  ################################################################################################################################


#' Performs fit step for a single SNP.
#'
#' Internal regression step for one SNP when learning U
#' Key options here include:
#'
#' carry_coeffs: this uses the previous iteration to initialize the current one. Thought it might provide a speed boost, but the difference seemed to be minimal.
#' glmnet: This uses the glmnet function; just exploring options
#' fastReg: This one was supposed to be faster. I don't know that it was at all; but not sure about serious results on it
#' ridge_L: this employs L2 as opposed to L1 regression  on L
#' fixed_ubiq: this allows for a ubiquitous factor by removing the L1 constraint from the first column.
#' @param x unscaled coefficients for single SNP across M studies
#' @param w Standard errors associated with row of x SNPs
#' @param W_c the decorrelation covariance matrix used in the GLS form of this software.
#' @param V Fixed factor matrix from previous iteration
#' @param option List of setting options (here, function accesses: traitSpecificVar, gls, carry_coeffs,regression_method,actively_calibrating_sparsity)
#' @param formerU Matrix of U from previous iteration to use in this one
#' @param r.v List of residual variances- for use only when accounting for study-specific variance
#'
#' @return a list containing: "l", the current learned row of U (k x 1) of weights for a given SNP, and the maximum possible sparsity for that run
#' @export
#'
#' @examples
#' source('../simulation/Generate_input.R');
#' data = generate_input(tau = tau);
#' F = data[[1]];
#' X = data[[3]];
#' W = data[[4]];
#'   option = list()
#'    option[['lambda1']] = 0.1;
#'   L_predict = fit_U(X, W, F, option);
#'
  one_fit <- function(x, w, W_c, V, option, formerU, r.v){
    if(option$traitSpecificVar && !is.null(r.v))
    {
      xp = x / sqrt(r.v);
      FactorMp = diag(1/sqrt(r.v)) %*% V;
    }else if(option$gls)
    {
      message("Don't use this option, we have a better way to do it...")
      #adjusting for covariance structure....
      #We need to make the covariance matrix is mostly empty...
      covar <- diag(1/w) %*% option$covar %*% diag(1/w) #rather than doing this for each SNP each time, we should just do it once for all SNPs
      #TODO- implement this later
      u_inv_t <- buildWhiteningMatrix(covar, blockify = FALSE) #already blockifyied the matrix...
      #this condition isn't necessarily reasonable, because adjustment might happen to force PD matrix.
      #might need to enforce this directly
      #stopifnot(all(diag(u_inv_t) == w)) #this isn't necessarily true. only if its a diagonal covar matri.
      xp <- adjustMatrixAsNeeded(x, covar, whitener=u_inv_t)
      FactorMp <- adjustMatrixAsNeeded(V, covar, whitener=u_inv_t)
    }
    else{
      xp = W_c %*% unlist(w * x);
      FactorMp = W_c %*% (diag(w) %*% V); #scaling this way is much smaller
    }
    #if we are doing the burn in....
    max.sparsity <- NULL
    ll = NULL
    penalty = NULL
    # Fit: xp' = FactorMp %*% l with l1 penalty on l -- |alpha1 * l|
    #Removed the transpose on X, getting the wrong dimensions
    dat_i = as.data.frame(cbind((xp), FactorMp));
    colnames(dat_i) = c('X', paste0('F', seq(1, ncol(FactorMp))));
    if(any(is.na(colnames(dat_i))) | any(is.na(paste0('F', seq(1,ncol(FactorMp))))))
    {
      message("In here, what is going on?")
      message("please debug.")
    }

    if(option[["carry_coeffs"]])
    {
      fit = penalized::penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=0,
                      positive = FALSE, standardize = option$std_coef, trace = FALSE, startbeta = formerU);
      l = penalized::coef(fit, 'all');
      ll = fit@loglik
      penalty = fit@penalty
    }else if (option[["regression_method"]] == "OLS") {
      fit <- RcppArmadillo::fastLmPure(as.matrix(FactorMp), xp)
      l = as.vector(coef(fit)) #check this, but I think its correct
      #ll = logLik(fit)
      ll = penLL(length(xp), xp - as.matrix(FactorMp) %*% fit$coefficients) # hack for now
    }
    else if(option[["regression_method"]] == "glmnet" & option[["fixed_ubiq"]])
      {
      penalties <- c(0, rep(1, (ncol(FactorMp) - 1)))
      fit <- glmnet::glmnet(x = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], y = xp, alpha = 1, lambda = option[['alpha1']],
                    intercept = FALSE, penalty.factor = penalties)
      #l = coef(fit, 'all')[,1][-1];
      l = glmnet::glmnetLASSO(dat_i, xp, ncol(FactorMp), option[['alpha1']], penalties)
      ll = loglik(fit)
    }
    else if(option[["regression_method"]] == "fastReg"){
      td <- t(dat_i[,paste0('F', seq(1,ncol(FactorMp)))])
      XTX <- td %*% dat_i[,paste0('F', seq(1,ncol(FactorMp)))]
      XTY <- td %*% xp
      fit <- elasticnet::elasticnet(XTX, XTY, lam2 = -1, lam1 = option[['alpha1']])
      ll = loglik(fit)
    }
    else if(option[["regression_method"]] == "ridge_L"){

      print("doing the ridge L")
      fit = penalized::penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 = 1e-15, lambda2=option[['alpha1']],
                      positive = FALSE, standardize = option$std_coef, trace = FALSE);
      l = penalized::coef(fit, 'all');
      ll = fit@loglik
      penalty = fit@penalty
    }
    else if( option[["regression_method"]] == "penalized" & option[["fixed_ubiq"]])
    {
      fit = penalized::penalized(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 =option[['alpha1']], lambda2=0,
                      positive = FALSE, standardize = option$std_coef, trace = FALSE) #maxiter = 30) #epsilon = option$epsilon)# #These limtations were hurting me...
      l = penalized::coef(fit, 'all')
      ll = fit@loglik
      penalty = fit@penalty
      #l = penalized::coef(fit, 'all')/penalized::weights(fit)
    }	else if(option$regression_method == "None"){
      l = c()

    }	else {
      fit = penalized::penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                      unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=0,
                      positive = FALSE, standardize = option$std_coef, trace = FALSE);
      l = penalized::coef(fit, 'all');
      ll = fit@loglik
      penalty = fit@penalty
    }
    if(option$actively_calibrating_sparsity) { max.sparsity <- rowiseMaxSparsity(xp, FactorMp)}
    return(list("l" = l, "max_sparse"=max.sparsity, "loglik" = ll, "penalty" = penalty))

  }


#how can I test and compare this?
  FitUGlobal <- function(X, W,W_c, V, option, formerU, r.v = NULL)
  {
    #it_U(X, W, W_c, initV, option)
    max.sparsity <- NA; penalty = NA; ll = NA; l = NA; penalty = NA
    long.x <- c(t(W_c) %*% t(X*W)) #stacks by SNP1, SNP2...
    #This multiplies each SNP row by the correction matrix
    weighted.copies <- lapply(1:nrow(X), function(i) t(W_c) %*% diag(W[i,]) %*% V)
    long.v <- Matrix::bdiag(weighted.copies) #weight this ish too you silly man.
    s = 1
    if(option$scale)
    {
      s <- getColScales(long.v)
      long.v <- unitScaleColumns(long.v, colnorms = s)
      #make into matrix for convenience:
      s = matrix(s, nrow = nrow(X), ncol = ncol(V),byrow = TRUE)
    }
    test.method = ""
    if(test.method == "preWeight")
    {
      s <- apply(V, 2, function(x) norm(x, "2"))
      s[s==0] <- 1 #so we don't divide by 0 on empty columns....
      scaled.v <- sweep(V,2,s,FUN="/")
      #stopifnot(norm(scaled.v[,1], "2") == 1)
      weighted.copies <- lapply(1:nrow(X), function(i) W_c %*% diag(W[i,]) %*% scaled.v)
      long.v <- Matrix::bdiag(weighted.copies) #weight this ish too you silly man.
    }
    if(test.method == "postWeight")
    {
      s = apply(long.v, 2, function(x) norm(x, "2"))
      s[s==0] <- 1 #so we don't divide by 0 on empty columns....
      long.v <- sweep(long.v,2,s,FUN="/")
      #long.v <- long.v / s
      #stopifnot(norm(long.v[,1], "2") == 1)
    }

    if(!is.null(formerU))
    {
      #CONCLLUSION: this does nothing to help the objective
      message("using old U")
      #formerU.scaled <- apply(formerU, 2, function(x) x/ norm(x, "2")) + 1e-15
      formerU.scaled <- formerU + 1e-15
      old.long.u <- c(t(formerU.scaled))
    }
    #formerU = NULL
    if(option$regression_method == "OLS")
    {
      fit <- RcppArmadillo::fastLmPure(as.matrix(long.v), long.x)
      l = as.vector(coef(fit)) #check this, but I think its correct
      #ll = logLik(fit)
      ll = penLL(length(long.x), long.x - as.matrix(long.v) %*% fit$coefficients) # hack for now
      l = matrix(l, nrow = nrow(X), ncol = ncol(V),byrow = TRUE)
    }else if( option$regression_method == "penalized")
    {
      if(is.null(formerU))
      {
        fit = penalized::penalized(response = long.x, penalized = as.matrix(long.v), lambda1 = option[['alpha1']],lambda2=0, unpenalized = ~0,
                                   positive = FALSE, standardize = option$std_coef, trace = FALSE)
      }else
      {
        fit = penalized::penalized(response = long.x, penalized = as.matrix(long.v), lambda1 = option[['alpha1']],lambda2=0, unpenalized = ~0,
                                   positive = FALSE, standardize = option$std_coef, trace = FALSE, startbeta = old.long.u)
      }

      l = penalized::coef(fit, 'all')
      ll = fit@loglik
      penalty = fit@penalty[1]
      l = matrix(l, nrow = nrow(X), ncol = ncol(V),byrow = TRUE)
    }else if(option$regression_method == "None"){
     l = c()
    }  else
    {
      message("Nothing implemented.")
    }

      #convert back to a U:

      if(option$actively_calibrating_sparsity) { max.sparsity <- rowiseMaxSparsity(as.matrix(long.x), as.matrix(long.v))}
      return(list("U" = l, "sparsity_space"=max.sparsity, "total.log.lik" = ll, "penalty" = penalty, "s"=s))
  }



  #U = fit_U(X, W, FactorM, option, ...);
  fit_U<- function(X, W,W_c, V, option, formerU, r.v = NULL){
    bad.cols <- c()
  	L = NULL
  	sparsity.est <- c()
  	running.log.lik <- 0
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
  		  l = one_fit(x, w, W_c, as.matrix(V), option, formerU[row,], r.v);
  		}		else {
  		  #l = one_fit(x, w, as.matrix(V), option, NULL, r.v);
  		  #tS = Sys.time();
  		  l = one_fit(x,w,W_c, as.matrix(V), option, NULL, r.v);
  		  #print(paste0('Updating Loading matrix takes ', round(difftime(Sys.time(), tS, units = "mins"), digits = 3)))
  		}
  		#Gettings some very bizarre errors about names and such, some weird hacks to work around it.
  		if(option$regression_method != "None")
  		{
    		if(length(l$l) == 1)
    		{
    		  l$l <- matrix(l$l)
    		  suppressMessages(names(l$l) <- "F1")
    		}
        if(is.null(L) | row == 1)
        {
            L = matrix(l$l, ncol = ncol(V))
        } else if(is.null(L) & is.null(l)){
            message("Matrix is empty...")
            return(list("U" = NULL, "sparsity_space"=sparsity.est))
        } else if(any(is.na(l$l)))
        {
          bad.cols <- c(bad.cols, which(is.na(l$l)))
          message("Dropping columns containing NA in burn-in, indicative of strong multi-colinearity")
          l$l[is.na(l$l)] <- 0
          L = suppressMessages(dplyr::bind_rows(L, l$l))
        }
        else{
          if(is.null(names(L)))
          {
            colnames(L) = paste0("F", 1:ncol(L))
          }
          #if(option$regression_method == "OLS")
         # {
            L = rbind(L, l$l)
          #}else
          #{
          #  L = suppressMessages(dplyr::bind_rows((L), (l$l)))
          #}
        }
    		}
    		if(option$actively_calibrating_sparsity) #we are in the process of burning in with OLS and calibrating sparsity
    		{
    		  sparsity.est <- c(sparsity.est, l$max_sparse)
    		}
    		running.log.lik <- running.log.lik + l$loglik
  	}
  	updateLog(paste0('Updating Loading matrix takes ', round(difftime(Sys.time(), tS, units = "mins"), digits = 3), ' min'), option);
    if(!is.null(L))
    {
      L <- as.matrix(L)
    }
  	return(list("U" = L, "sparsity_space"=sparsity.est, "redundant_cols" = unique(bad.cols), "total.log.lik"=running.log.lik))
  }


#Should be working fine... need to do some testing though!
  fitUParallel <- function(X, W, W_c, V, option, formerU){

    L = NULL
    tS = Sys.time()
    #cl <- parallel::makeCluster(option[["ncores"]])#, outfile = "TEST.txt")
    iterations <- nrow(X)
    #doParallel::registerDoParallel(cl)
    doParallel:registerDoParallel(parallel::makeCluster(option[["ncores"]]))
    #6/13/2022
    #We split the tasks across cores, so multiple regression steps per core; you aren't switching nearly as much.
    n.per.group <- ceiling(nrow(X)/option$nsplits)
    oo <- n.per.group - 1
    split_lines <- lapply(1:(option$nsplits-1), function(x) (x*n.per.group - oo) :(x*n.per.group))
    split_lines[[option$nsplits]] <- (split_lines[[(option$nsplits-1)]][n.per.group]+1):nrow(X)
    L <- foreach(rows = split_lines, .combine = 'dplyr::bind_rows', .packages = c('penalized', 'dplyr')) %dopar% {
      sub_l <- NULL
      for(row in rows)
      {
        x = X[row, ];
        w = W[row, ];

        xp = W_c %*% matrix(w * x, nrow = nrow(V), ncol = 1); #elementwise multiplication
        FactorMp = W_c %*% diag(w) %*% v;  #same as below
        dat_i = data.frame(cbind(xp, FactorMp));
        colnames(dat_i) = c('X', paste0('F', seq(1, ncol(FactorMp))));
        fit = penalized::penalized(response = X, penalized = dat_i[,paste0('F', seq(1,ncol(FactorMp)))], data=dat_i,
                        unpenalized = ~0, lambda1 = option[['alpha1']], lambda2=0,
                        positive = FALSE, standardize = option$std_coef, trace = FALSE,epsilon = option$epsilon);
        sub_l <- dplyr::bind_rows(sub_l, penalized::coef(fit, 'all'))
      }
      sub_l
    }
    doParallel::stopImplicitCluster()
   #parallel::stopCluster(cl)
    updateLog(paste0('Updating Loading matrix takes ',  as.character(round(difftime(Sys.time(), tS, units = "mins"), digits = 3), ' min')), option);
    return(as.matrix(L))
  }



#' Wrapper for the U-fitting step
#'
#' @param X the matrix to be decomposed (N x T, N is the number of data points, T is the number of features)
#' @param W the weight matrix, same size as in X, containing standard errors of the GWAS estimates
#' @param W_c the whitening matrix, associated with covariance structure, a M x M
#' @param V learned factor matrix  (M x K, M is the number of features, K is the number of factors)
#' @param option list with parameters including lambda1, the l1 penalty parameter for the U matrix
#' (relevant arguments include: ncores,)
#' @param ... arguments passed on to FitU
#'
#' @return A list containing the return matrix U, the distribution of maximum sparsity parameters, and any columns dropped because they contained NAs
#' @export
#'
  FitUWrapper <- function(X,W,W_c,V, option, prevU = NULL,...)
  {
    if(option[['ncores']] > 1) #This is not working at all. Can't tell you why. But its not. Need to spend some time debugging at some point.
    {
      log_print("Fitting L in parallel")
      U = fitUParallel(X, W,W_c, V, option, formerU = prevU); #preL is by default Null, unless yo specify!

    }else
    {
      #U = fit_U(X, W, W_c, V, option,formerU = prevU, ...);
      U = FitUGlobal(X, W,W_c, V, option,formerU = prevU, ...);

    }
    return(U)
  }


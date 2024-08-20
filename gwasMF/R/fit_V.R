################################################################################################################################
#Modified 4/26 by Ashton Omdahl, based on previous work by Yuan He
##
## Input:
##		- X:
##		- L:
##		- W:
##		- option: a list with parameters including lambda1, the l1 penalty parameter for the factor matrix (F)
##
## Return:
##		- A non-negative factor matrix that minimize_F ||(X - LF') .* W||_F^2 + lambda1*|F|_1, F is non-negative
##
## Example of usage:
##
## source('../simulation/Generate_input.R');
## data = generate_input(tau = tau);
## L = data[[2]];
## X = data[[3]];
## W = data[[4]];
## option = list();
## option[['lambda1']] = 0.1;
## F_predict = fit_F(X, W, L, option);
##
################################################################################################################################


#' Update the factor matrix V during the alternating least squares
#'
#' @param X the matrix to be decomposed (N x T, N is the number of data points, T is the number of features)
#' @param W the weight matrix, same size as in X
#' @param U learned loading matrix  (N x K, N is the number of data points, K is the number of factors)
#' @param option list of all option settings (accesses V, LD, gls,
#' actively_calibrating_sparsity, posF, regression_method, epsilon,traitSpecificVar,intercept_ubiq)
#' @param formerV Optional argument to pass in previous iteration of V as prior to next one.
#'
#' @return List containing the new V ("V") and the sparsity upper limits for each regression steps ("sparsity_space")
#' @export
#'
fit_V <- function(X, W, U, option, formerV = NULL){
	tStart   = Sys.time();
	FactorM  = NULL;
	r.v <- c()
	lambda_list <- c()
	sparsity.est <- c()
	running.loglik <- 0
	## fit each factor one by one -- because of the element-wise multiplication from weights!
	for(col in seq(1, ncol(X))){
	  ll = NULL
		if(option$V > 0)
		{
			svMisc::progress(col, ncol(X), progress.bar = TRUE)
		}
    	x = X[, col];
	    w = W[, col];
        ## weight the equations on both sides
    if(option$gls)
    {

      covar <- diag(1/w) %*% option$LD %*% diag(1/w) #rather than doing this for each SNP each time, we should just do it once for all SNPs
      u_inv_t <- buildWhiteningMatrix(covar, blockify = FALSE) #already blockifyied the matrix...
      xp <- adjustMatrixAsNeeded(x, covar, whitener=u_inv_t)
      Lp <- adjustMatrixAsNeeded(U, covar, whitener=u_inv_t)
    }else
    {
      xp = unlist(w) * x;
      Lp = as.matrix(unlist(w)* U); #CONFIRMEd: this has desired behavior.
    }

    if(all(Lp == 0))
    {
      print("all 0 case, beware.")
    }
    if(option$actively_calibrating_sparsity) { sparsity.est <- c(sparsity.est, rowiseMaxSparsity(xp, Lp, fixed_first = option$fixed_ubiq))}
		## Fit: xp = U %*% f with l1 penalty on f -- |lambda1 * f|
		dat_i = as.data.frame(cbind(xp, Lp));
		colnames(dat_i) = c('X', paste0('F', seq(1, ncol(Lp))));
		if(any(is.na(colnames(dat_i))) | any(is.na(paste0('F', seq(2,ncol(Lp))))))
		{
		  message("In here, what is going on?")
		  message("please debug.")
		}
		# unpenalized = ~0 - suppress the intercept
		# positive = TRUE  - constrains the estimated regression coefficients of all penalized covariates to be non-negative
		# lambda2=1e-10    - avoid singularity in fitting
		if(option[["carry_coeffs"]])
		{
		  fit = penalized::penalized(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(Lp)))], data=dat_i,
		                  unpenalized = ~0, lambda1 = option[['lambda1']], lambda2=0,
		                  positive = option$posF, standardize = option$std_coef, trace = FALSE, startbeta = formerV[col,], maxiter=50 )
		  f = penalized::coef(fit, 'all')
		  ll = fit@loglik
		} else if(option[["regression_method"]] == "OLS")
		{
		  fit <- lm(xp~Lp + 0)
		  f = coef(fit) #check this, but I think its correct
		}
		else if(option[["regression_method"]] == "glmnet" & option[["fixed_ubiq"]])
    {
				penalties <- c(0, rep(1, (ncol(L) - 1)))
				f = glmnetLASSO(dat_i, xp, ncol(U), option[['lambda1']], penalties)
   	} else if(option[["fixed_ubiq"]] & option[["regression_method"]] == "penalized")
		{
		  #lambdas <- c(0, rep(option[['lambda1']], (ncol(Lp) - 1)))
   	  lambdas <- c(rep(option[['lambda1']], (ncol(Lp) - 1)))#no lasso on that first column
		  #lambdas.2 <- c(0.5, rep(0, (ncol(Lp) - 1)))
   	  if(length(lambdas) == 0)
   	  {
   	    fit <- RcppArmadillo::fastLmPure(as.matrix(Lp), xp)
   	    f = as.vector(coef(fit)) #check this, but I think its correct
   	    #ll = logLik(fit)
   	    ll = penLL(length(xp), xp - as.matrix(Lp) %*% fit$coefficients) # hack for now
   	  }else
   	  {
   	    fit <- penalized::penalized(response = X, penalized = dat_i[,paste0('F', seq(2, ncol(Lp)))], data=dat_i,
   	                                unpenalized = ~F1 + 0, lambda1 =option[['lambda1']], lambda2=0,
   	                                positive = option$posF, standardize = option$std_coef, trace = FALSE) #, epsilon = option$epsilon, maxiter= 50)


   	    f = penalized::coef(fit, 'all')
   	    ll = fit@loglik
   	  }
		  #f= penalized::coef(fit, 'all')/penalized::weights(fit)
		} else if(option$regression_method == "None")
		{
		  f = c()
		}
		else {
		  fit = penalized::penalized(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(Lp)))], data=dat_i,
		                  unpenalized = ~0, lambda1 =option[['lambda1']], lambda2=0,
		                  positive = option$posF, standardize = option$std_coef, trace = FALSE, maxiter= 50)
		  f = penalized::coef(fit, 'all')
		  ll = fit@loglik
		}

		if(option$traitSpecificVar)
		{
		  n = nrow(dat_i[,paste0('F', seq(1, ncol(Lp)))])
		  p = ncol(dat_i[,paste0('F', seq(1, ncol(Lp)))])
		  res.var <- (1/(n-p)) * fit@residuals %*% fit@residuals
		  r.v <- c(r.v, res.var)
		}
		FactorM = rbind(FactorM, f);
		running.loglik <- running.loglik + ll

		#track the objective, if set

	}

	updateLog(paste0('Updating Factor matrix takes ',  round(difftime(Sys.time(), tStart, units = "mins"), digits = 3), ' min'), option)
	if(option$traitSpecificVar)
	{
	  ret <- list()
	  ret$r.v <- r.v
	  ret$mat <- FactorM
	  ret$sparsity_space = sparsity.est
	  return(ret)
	}


  #  if(option$intercept_ubiq) #Force the ubiquitous term to be fixed. Run with fixed_ubiq? doesn't really matter{
  #  {
  #    message("Not sure why'd you do this, but okay...")
  #    cor_struct <- cor(X)
  #    svd <- svd(cor_struct, nu = 1) #fortunately its symmetric, so  U and V are the same here!
  #    ones <- sign(svd$u)
  #      FactorM[,1] <- ones
  #  }

    return(list("V" = FactorM, "sparsity_space"=sparsity.est, "total.log.lik" = running.loglik))
}

FitVGlobal <- function(X, W, W_c, U, option, formerV = NULL, reg.elements=NULL)
{
  #avail.cores <- parallel::detectCores() -1
  K <- ncol(U)
  M = ncol(X)
  N=nrow(X)
   sparsity.est <- NA; penalty = NA; ll = NA; l = NA; FactorM = c(); penalty = NA
   if(is.null(reg.elements))
   {
     long.x <- c(t(W_c) %*% t(X*W)) #Stacked by SNP
     if(option$std_y) { long.x <- mleStdVector(long.x)}
     #Weights to multiple the reordered elements of U by
     joined.weights <- lapply(1:nrow(X), function(i) t(W_c) %*% (diag(W[i,])))
     #print(paste0("made joined weight: ", pryr::mem_used()))
   }else
   {
     long.x <- reg.elements$long.x
     joined.weights <-  reg.elements$joined.weights
   }
   s = 1
   if(FALSE) #Old way of building U
   {
     #Old way- new way is faster
     nsnps = min(1000, nrow(U))
     interval <- floor(nrow(X)/nsnps)
     extra.count <- nrow(X) - interval * nsnps
     blocks <- lapply(1:interval, function(i) ((i-1)*nsnps + 1):(i*nsnps))
     if(extra.count != 0)
     {
       blocks[[length(blocks) + 1]] <- (interval * nsnps + 1) : (interval * nsnps + extra.count)
     }


     print("warning-old way....")
     tictoc::tic()
     all.pieces <- lapply(blocks, function(x) stackUs(x, M,K, joined.weights, U, norm = option$scale))
     #all.pieces <- parallel:mclapply(blocks, function(x) stackUs(x, M,K, joined.weights, U, norm = option$scale), mc.cores=avail.cores)
     long.u <- do.call("rbind",all.pieces )
     # long.u <- do.call("rbind",lapply(blocks, function(x) stackUs(x, M,K, joined.weights, U, norm = option$scale)) )
     tictoc::toc()
     if(option$scale)
     {
       #s.tmp <- FrobScale(long.u)
       #s <- s.tmp$s; long.u <- s.tmp$m.scaled
       #add across all blocks
       all.pieces <- lapply(blocks, function(x) stackUs(x, M,K, joined.weights, U, norm = option$scale))
       long.u <- do.call("rbind", lapply(all.pieces, function(x) x$matrices))
       s <- max(sapply(1:M, function(i) sum(sapply(all.pieces, function(x) x$norm.list[[i]]))))
       #s <- max(lapply(all.pieces, function(x) x$norm.list))
       long.u <- long.u/s
     }

   }

  long.u <- completeUExpand(joined.weights, nrow(U),M,K,U)
  #message("New u added: ", lobstr::mem_used())
  nopen.cols <- sapply(1:ncol(long.u), function(x) x %% K)

  #special case this logic failes on
  if(K == 1 & option$fixed_ubiq & option$regression_method == "penalized")
  {
    #we are just down to the last column
    option$regression_method = "OLS"
    #Just bypass this.
  }

  dat_i = long.u
  unpen <- dat_i[,nopen.cols == 1]
  pen <- dat_i #include everythign if all penalized
 if(option$fixed_ubiq & option$regression_method == "penalized")
{
   if(option$fixed_ubiq)
   {
     #only penalize the correct columns.
     pen <- dat_i[,nopen.cols != 1]
   }

   if(!is.null(formerV))
   {
    #old.long.v <- c(t(formerV))
     fixed.starts <- formerV[,1]
     message("doing it this way..")
     penalized.starts <-  c(t(formerV[,-1]))
     fit <- penalized::penalized(response = as.matrix(long.x), penalized =as.matrix(pen),
                                 unpenalized = ~as.matrix(unpen) + 0, lambda1 =option[['lambda1']], lambda2=0,
                                 positive = option$posF, standardize = option$std_coef, trace = FALSE,startbeta = penalized.starts, startgamma = fixed.starts ) #, epsilon = option$epsilon, maxiter= 50)

   }else
   {
     fit <- penalized::penalized(response = as.matrix(long.x), penalized =as.matrix(pen),
                                 unpenalized = ~as.matrix(unpen) + 0, lambda1 =option[['lambda1']], lambda2=0,
                                 positive = option$posF, standardize = option$std_coef, trace = FALSE) #, epsilon = option$epsilon, maxiter= 50)
   }
   f = penalized::coef(fit, 'all')
      ll = fit@loglik
      #F 1-M are the 1st factor
      FactorM = cbind(f[1:M],matrix(f[(M+1):(M*K)], nrow = M,byrow = TRUE))
      penalty = fit@penalty[1]
} else if( option$regression_method == "glmnet")
{
  lasso.active.set <- rep(1, ncol(long.u))
  if( option$fixed_ubiq) {lasso.active.set[nopen.cols == 1] <- 0 }

 if(is.na(option$lambda1)) #we are still parameter searching
 {
   fit <- glmnet::glmnet(x = long.u, y = long.x, family = "gaussian", alpha = 1,
                        intercept = FALSE, standardize = option$std_coef,
                        penalty.factor = lasso.active.set,  trace.it=1) #lambda = option[['alpha1']],
   #use the pryr:mem_used() and mem_changed() functions to track this.
   lambda.list <- fit$lambda
   penalty <- fit$penalty
   if(FALSE)
   {   bic.list <- c()
     for(i in 1:length(lambda.list))
     {
       #convert u into a matrix we can work with
       #have to remove the intercept term.
       v.curr = matrix(coef(fit, s = lambda.list[i])[-1], nrow = ncol(X), ncol = ncol(U),byrow = TRUE)
       #calculate a BIC for each setting
       bic.list <- c(bic.list, CalcMatrixBIC.loglikGLOBALversion(X,W,U,v.curr, W_cov = W_c, which.learning = "V", df = fit$df[i],fixed_first=option$fixed_ubiq))

     }
   }
   #message("Estimating BIC... ", pryr::mem_used())
   bic.list <- calculateBIC(fit, long.u, long.x, option$bic.var)

  #print(paste0("Calculated BIC: ", pryr::mem_used()))
   #print(pryr::mem_used())
   #bic.list <- BICglm(fit, option$bic.var)
   #bic.list <- sklearnBIC(fit,long.u,long.x, bic.mode =  option$bic.var)
   #bic.list <- ZouBIC(fit, long.u, long.x)
   #pick the best one
   min.index <- which.min(bic.list)
   v.ret = matrix(coef(fit, s = lambda.list[min.index])[-1], nrow = ncol(X), ncol = ncol(U),byrow = TRUE)

   stopifnot(!is.unsorted(-fit$lambda)) #make sure its ordered
   #grab the previous 25 and the next 25
   upper.next <- NA; lower.next <- NA
   if(min.index > 25)
   {
     upper.next <- fit$lambda[min.index - 25]
   }
   if((min.index + 25) < length(fit$lambda))
   {
     lower.next <- fit$lambda[min.index + 25]
   }
   return(list("V" =  pm(v.ret),"lambda.sel"=fit$lambda[min.index],"bic"= bic.list[min.index], "sparsity_space"=max(fit$lambda),
               "total.log.lik" = NA, "penalty" = penalty, "s"=s, "next.upper" = upper.next, "next.lower" = lower.next))

 }else
 {
   fit <- glmnet::glmnet(x = long.u, y = long.x, family = "gaussian", alpha = 1,
                        intercept = FALSE, standardize = option$std_coef, penalty.factor = lasso.active.set,
                        lambda=option$lambda1, trace.it=1)
   #print(paste0("Fitting done- every first entry is 0: ", pryr::mem_used()))
   penalty <- fit$penalty
   v.curr = matrix(coef(fit, s = option$lambda1)[-1], nrow = ncol(X), ncol = ncol(U),byrow = TRUE)
   return(list("V" = pm(v.curr), "sparsity_space"=max(fit$lambda), "total.log.lik" = NA, "penalty" = NA, "s"=s,"SSE"=deviance(fit)))
 }

  #ProposeNewSparsityParams(bic.list, fit$lambda, (fit$lambda), 2, n.points = 20) #need to modify this to allow for conditions.....
  #Check: is in the realm of those picked by CV?


} else if(option$regression_method == "OLS")
  {
    fit <- RcppArmadillo::fastLmPure(as.matrix(long.u), long.x)
    f = as.vector(coef(fit)) #check this, but I think its correct
    FactorM = matrix(f, nrow = M,byrow = TRUE)
    ll = penLL(length(long.x), long.x - as.matrix(long.u) %*% fit$coefficients) # hack for now)
    #here we don't even pass in the 1st factor, so not an issue.
}	else if(option$regression_method == "None"){
  FactorM = c()
  }  else
    {
      message("Not implemented. Crash and burn.")
    }
  if(option$actively_calibrating_sparsity) { sparsity.est <- rowiseMaxSparsity(long.x, pen, fixed_first = FALSE)}

  return(list("V" = pm(FactorM), "sparsity_space"=sparsity.est, "total.log.lik" = ll, "penalty" = penalty, "s"=s))

}
#' Wrapper for the V function, not serving a purpose at the moment.
#'
#' @param X the matrix to be decomposed (N x T, N is the number of data points, T is the number of features)
#' @param W the weight matrix, same size as in X, containing standard errors of the GWAS estimates
#' @param U learned loading matrix  (N x K, N is the number of data points, K is the number of factors)
#' @param option list of all option settings (accesses V, LD, gls,
#' actively_calibrating_sparsity, posF, regression_method, epsilon,traitSpecificVar,intercept_ubiq)
#' @param formerV Optional argument to pass in previous iteration of V as prior to next one.
#'
#' @return List containing the new V ("V") and the sparsity upper limits for each regression steps ("sparsity_space")
#' @export
#'
FitVWrapper <- function(X, W,W_c, U, option, formerV = NULL,...)
{
  FitVGlobal(X, W, W_c, as.matrix(U), option, formerV = formerV,...)
}




##HELPER FUNCTION for quickly stacking vs sparsely

#' Stack weighted elements of U to learn V. Could be helpful for paralellizing
#'
#' @param nsnps the chunk size to analyze
#' @param M the number of studies
#' @param K current K
#' @param joined.weights a list of the covariance matrix times the diagonal weightings matrices
#' @param U current u
#'
#' @return an nm x mk sparse matrix
#' @export
stackUs <- function(nsnps, M, K, joined.weights, U, norm = FALSE)
{
  if(norm){norm.list <- rep(0, M)}
  test <- do.call("rbind", lapply(0:(length(nsnps)-1), function(n) matrix(c(sapply(1:M, function(x) rep((x + n*M),K)),
                                                                            1:(M * K),
                                                                            rep(U[(n + min(nsnps)),],M)),
                                                                          ncol = 3)))



  ret.dat.alt <- Matrix::bdiag(joined.weights[nsnps]) %*% Matrix::sparseMatrix(i = test[,1], j = test[,2], x = test[,3], dims = c(length(nsnps) * M, M*K))
  if(norm)
  {
    return(list("matrices"=ret.dat.alt, "norm.list"=lapply(norm.list, sqrt)))
  }
  return(ret.dat.alt)
}

#An alternative to compare
completeUExpand <- function(joined.weights, N,M,K,U)
{
  test <- do.call("rbind", lapply(1:N, function(n) matrix(c(sapply(1:M, function(m) rep(m+(n-1)*M,K)),
                                                            1:(M * K),
                                                            rep(U[n,],M)),
                                                          ncol = 3)))
  ret.dat.alt <- Matrix::sparseMatrix(i = test[,1], j = test[,2], x = test[,3], dims = c(N * M, M*K))
  Matrix::bdiag(joined.weights) %*% ret.dat.alt
}




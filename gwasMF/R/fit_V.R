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
	obj.updates <- c() #tracking the objective at each step.
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
      #TODO: simplify. Just need C, not both S and C here.
      #adjusting for covariance structure.... here it is LD in U
      #We need to make the covariance matrix is mostly empty...
      covar <- diag(1/w) %*% option$LD %*% diag(1/w) #rather than doing this for each SNP each time, we should just do it once for all SNPs
      #TODO- implement this later
      u_inv_t <- buildWhiteningMatrix(covar, blockify = FALSE) #already blockifyied the matrix...
      #this condition isn't necessarily reasonable, because adjustment might happen to force PD matrix.
      #might need to enforce this directly
      #stopifnot(all(diag(u_inv_t) == w)) #this isn't necessarily true. only if its a diagonal covar matri.
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
   	                                unpenalized = ~F1 + 0, lambda1 =option[['lambda1']], lambda2=1e-10,
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
#if(option$debug) {
#  obj.updates <- c(obj.updates, compute_obj(X,W,W_c,U,))
#  (X, W, U, option, formerV = NULL)
#  compute_obj <- function(X, W, W_c, L, FactorM, option, decomp = FALSE, loglik = TRUE, globalLL=FALSE)
#}

  #  if(option$intercept_ubiq) #Force the ubiquitous term to be fixed. Run with fixed_ubiq? doesn't really matter{
  #  {
  #    message("Not sure why'd you do this, but okay...")
  #    cor_struct <- cor(X)
  #    svd <- svd(cor_struct, nu = 1) #fortunately its symmetric, so  U and V are the same here!
  #    ones <- sign(svd$u)
  #      FactorM[,1] <- ones
  #  }

    return(list("V" = FactorM, "sparsity_space"=sparsity.est, "total.log.lik" = running.loglik, "obj.changes" = obj.updates))
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
FitVWrapper <- function(X, W, U, option, formerV = NULL)
{
  #HERE?
  fit_V(X, W, as.matrix(U), option, formerV = formerV)
}

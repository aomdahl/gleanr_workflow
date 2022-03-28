################################################################################################################################
#Modified 4/26 by Ashton Omdahl
## Update the factor matrix F during the alternative optimization
##
## Input:
##		- X: the matrix to be decomposed (N x T, N is the number of data points, T is the number of features)
##		- L: learned loading matrix  (N x K, N is the number of data points, K is the number of factors)
##		- W: the weight matrix, same size as in X
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

suppressWarnings(library(penalized))

fit_F <- function(X, W, L, option, formerF = NULL){
	tStart   = Sys.time();
	FactorM  = NULL;
	r.v <- c()
	lambda_list <- c()
	## fit each factor one by one -- because of the element-wise multiplication from weights!
	for(col in seq(1, ncol(X))){
        
      #x = X[, col, with = FALSE];
      #w = W[, col, with = FALSE];
    x = X[, col];
	    w = W[, col];
        ## weight the equations on both sides
		xp = unlist(w) * x;
		Lp = unlist(w)* L;

		## Fit: xp = L %*% f with l1 penalty on f -- |lambda1 * f|
		dat_i = as.data.frame(cbind(xp, Lp));
		colnames(dat_i) = c('X', paste0('F', seq(1, ncol(Lp))));

		# unpenalized = ~0 - suppress the intercept
		# positive = TRUE  - constrains the estimated regression coefficients of all penalized covariates to be non-negative
		# lambda2=1e-10    - avoid singularity in fitting
		if(option$calibrate_sparsity){
		  f <- recommendRange(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(Lp)))], data=dat_i,
		                      unpenalized = ~0, lambda1 =option[['lambda1']], lambda2=1e-10,
		                      positive = option$posF)
		  
		}	else if(option[["reweighted"]])
		{
		  fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(Lp)))], data=dat_i,
		                  unpenalized = ~0, lambda1 = option[['lambda1']], lambda2=1e-10,
		                  positive = option$posF, standardize = FALSE, trace = FALSE, startbeta = formerF[col,] )
		  f = coef(fit, 'all')
		  #fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(2, ncol(Lp)))], data=dat_i,
		                  #unpenalized = ~0 + dat_i[,1], lambda1 = option[['lambda1']], lambda2=1e-10,
		                  #positive = FALSE, standardize = FALSE, trace = FALSE, startbeta = formerF[col,] )
		} else if(option[["glmnet"]])
		{
		  glm_fit <- glmnet(x = dat_i[,paste0('F', seq(1, ncol(Lp)))], y = xp, alpha = 1, lambda = c(option[['lambda1']]), intercept = FALSE)
		  f = coef(glm_fit, 'all')[,1][-1]
		} else if(option[["fixed_ubiq"]])
		{
		  lambdas <- c(0, rep(option[['lambda1']], (ncol(Lp) - 1))) #no lasso on that first column
		  fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(Lp)))], data=dat_i,
		                  unpenalized = ~0, lambda1 =lambdas, lambda2=1e-10,
		                  positive = option$posF, standardize = FALSE, trace = FALSE)
		  #glm_fit <- glmnet(x = dat_i[,paste0('F', seq(1, ncol(Lp)))], y = xp, alpha = 1, lambda = lambdas, intercept = FALSE)
		  
		  f = coef(fit, 'all')
		} else {
		  fit = penalized(response = X, penalized = dat_i[,paste0('F', seq(1, ncol(Lp)))], data=dat_i,
		                  unpenalized = ~0, lambda1 =option[['lambda1']], lambda2=1e-10,
		                  positive = option$posF, standardize = FALSE, trace = FALSE)
		  f = coef(fit, 'all')
		}
		
		if(option$traitSpecificVar)
		{
		  n = nrow(dat_i[,paste0('F', seq(1, ncol(Lp)))])
		  p = ncol(dat_i[,paste0('F', seq(1, ncol(Lp)))])
		  res.var <- (1/(n-p)) * fit@residuals %*% fit@residuals
		  r.v <- c(r.v, res.var)
		}
		FactorM = rbind(FactorM, f);
		#print(col)
	}

	tEnd = Sys.time();
	message(paste0('Updating Factor matrix takes ', round((tEnd - tStart)/60, 2), 'min'));
	if(option$traitSpecificVar)
	{
	  ret <- list()
	  ret$r.v <- r.v
	  ret$mat <- FactorM
	  return(ret)
	}

    if(option$intercept_ubiq) #Force the ubiquitous term to be fixed. Run with fixed_ubiq? doesn't really matter{
    {
        cor_struct <- cor(X)
      svd <- svd(cor_struct, nu = 1) #fortunately its symmetric, so  U and V are the same here!
      ones <- sign(svd$u) 
        FactorM[,1] <- ones
	}
    return(FactorM)
}

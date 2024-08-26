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
                    intercept = FALSE, penalty.factor = penalties, trace.it=1)
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
  FitUGlobal <- function(X, W,W_c, V, option, formerU, r.v = NULL,reg.elements=NULL)
  {
    #it_U(X, W, W_c, initV, option)
    #print(paste0("Starting global U fit: ", pryr::mem_used()))
    #TweaK this:
    #glmnet::glmnet.control(fdev = 1e-10)

    max.sparsity <- NA; penalty = NA; ll = NA; l = NA; sse=NA
    if(is.null(reg.elements))
    {
      long.x <- c(t(W_c) %*% t(X*W)) #stacks by SNP1, SNP2... update so not needing to be passed each time.
      joined.weights <- lapply(1:nrow(X), function(i) Matrix::Matrix(t(W_c) %*% (diag(W[i,]))))
      if(option$std_y) { long.x <- mleStdVector(long.x)}
    }else
    {
      long.x <- reg.elements$long.x
      joined.weights <- reg.elements$joined.weights
    }
    if(option$std_y) { long.x <- mleStdVector(long.x)}

    #This multiplies each SNP row by the correction matrix
    #weighted.copies <- lapply(1:nrow(X), function(i) t(W_c) %*% diag(W[i,]) %*% V)
    #long.v <- Matrix::bdiag(weighted.copies) #weight this ish too you silly man.
    #TODO- use the data in reg.elements so we aren't repeating work here- multiply all sparse elemnts of regg.elements$mat by V
    #fastest is 3rd one. ya boi.
    #long.v <- Matrix::bdiag(lapply(1:nrow(X), function(i) t(W_c) %*% diag(W[i,]) %*% V)) #this is quite fast actually
    #long.v <- Matrix::bdiag(reg.elements$joined.weights) %*% Matrix::bdiag(replicate(nrow(X),V, simplify=FALSE))
    long.v <- Matrix::bdiag(lapply(joined.weights, function(x) x %*% V))
    #see which of the 3 above is faster and more memory efficient for use going forward.
    #print(paste0("Built sparse matrix: ", pryr::mem_used()))
    #message("This matrix object takes up: ", object.size(long.v))
  s = 1
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
    }else if( option$regression_method == "glmnet")
    {
      if(is.na(option$alpha1)) #we are searching for parameters, phase 1 of the work.
      {
        fit <- glmnet::glmnet(x = long.v, y = long.x, family = "gaussian", alpha = 1,
                             intercept = FALSE, standardize = option$std_coef,nlambda = 100, trace.it=1) #lambda = option[['alpha1']],
        #print(paste0("Just fit the regression: ", pryr::mem_used()))
        penalty <- fit$penalty
        alpha.list <- fit$lambda
        #message("SAVING, and Using BIC to select out the best one...")
        #save(X,W,W_cfit, long.v, long.x, file="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/DEBUG.bic.RData" )
      if(FALSE)
      {
        all.preds <- predict(fit, newx = long.v)
        #Is all preds closer to the scaled or unscaled version of longx?
        #stopifnot(sum((all.preds[,100] - long.x)^2) < sum((all.preds[,100] - mleStdVector(long.x))^2))

        #Documentation says: " standardizes y to have unit variance (using 1/n rather than 1/(n-1) formula)
        #before computing its lambda sequence (and then unstandardizes the resulting coefficients)"
        #So as we see here, use the UNSTANDARDIZED version.

        #Almost CERTAINLY the unscaled version The distance on the scaled is much bigger
        bic.list.complete<- list()
        bic.list.complete[["bic.std"]] <- c();bic.list.complete[["bic.std_prev"]] <- c()
        bic.list.complete[["bic.std.scaled"]] <- c(); # not doing this b/c same as prev, too much work.bic.list.complete[["bic.std_prev.scaled"]] <- c()
        bic.list.complete[["aic"]] <- c();bic.list.complete[["aic.scaled"]] <- c();
        bic.list.complete[["avg"]] <- c(); bic.list.complete[["avg.scaled"]] <- c()
        bic.list.complete[["zou"]] <- c(); bic.list.complete[["zou.scaled"]] <- c()
        bic.list.complete[["dev"]] <- c(); bic.list.complete[["zou.unbiased"]] <- c();
          for(i in 1:length(alpha.list))
          {

            #May 15 updaed
            bic.list.complete[["bic.std"]] <- c(bic.list.complete[["bic.std"]], (-2*penLLSimp(nrow(X), ncol(X), long.x - all.preds[,i]) + log(length(long.x)) * fit$df[i]))
            bic.list.complete[["aic"]] <- c(bic.list.complete[["aic"]], (-2*penLLSimp(nrow(X), ncol(X), long.x - all.preds[,i]) + 2 * fit$df[i]))
            bic.list.complete[["avg"]] <- c(bic.list.complete[["avg"]], (-2*penLLSimp(nrow(X), ncol(X), long.x - all.preds[,i]) + mean(c(2,log(length(long.x)))) * fit$df[i]))
            #convert u into a matrix we can work with
            #have to remove the intercept term.
            u.curr = matrix(coef(fit, s = alpha.list[i])[-1], nrow = nrow(X), ncol = ncol(V),byrow = TRUE)
            bic.list.complete[["bic.std_prev"]] <- c( bic.list.complete[["bic.std_prev"]], CalcMatrixBIC.loglikGLOBALversion(X,W,u.curr,V, W_cov = W_c,
                                                                      which.learning = "U", df = fit$df[i],fixed_first=FALSE, scale.explanatory = FALSE))

            #Now the scaled version of each, which is what I think it really should be......
            scaled.long.x <-  mleStdVector(long.x)
            bic.list.complete[["bic.std.scaled"]] <- c(bic.list.complete[["bic.std.scaled"]], (-2*penLLSimp(nrow(X), ncol(X), scaled.long.x - all.preds[,i]) + log(length(scaled.long.x)) * fit$df[i]))
            bic.list.complete[["aic.scaled"]] <- c( bic.list.complete[["aic.scaled"]], (-2*penLLSimp(nrow(X), ncol(X), scaled.long.x - all.preds[,i]) + 2 * fit$df[i]))
            bic.list.complete[["avg.scaled"]] <- c(bic.list.complete[["avg.scaled"]], (-2*penLLSimp(nrow(X), ncol(X), scaled.long.x - all.preds[,i]) + mean(c(2,log(length(scaled.long.x)))) * fit$df[i]))
               #Verified- residuals are the SAME!
            #resid.new <- long.x - all.preds[,i]
            #resid.old <- calcGlobalResiduals(X, W, u.curr, V, W_cov = W_c)

          }
        bic.list.complete[["zou"]] <- ZouBIC(fit, long.v, long.x)
        bic.list.complete[["zou.unbiased"]] <- ZouBIC(fit, long.v, long.x,  var.meth = "unbiased")
        bic.list.complete[["dev"]] <- BICglm(fit)
        bic.list.complete[["dev.e"]] <- BICglm(fit,bic.mode = "ebic")
      }

        #bic.list <- BICglm(fit, option$bic.var)
        #message("Calc BIC")
        #print(dim(long.v))
        #print(length(long.x))
        bic.list <- calculateBIC(fit, long.v, long.x, option$bic.var)
        #print(paste0("Calculated BIC: ", pryr::mem_used()))
        #print(pryr::mem_used())
        #sklearn extended was default previously?
        #bic.list <- sklearnBIC(fit,long.v,long.x, bic.mode =  option$bic.var)
        #save(bic.list, fit, long.x, long.v,file="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/DEBUG.bic.RData" )
        #bic.list <- ZouBIC(fit, long.v, long.x) #trying this...



      if(FALSE)
      {
        #Trying generalized hbic for high dimensional case, Wang et al.
        c = log(log(length(long.x)))
        n <- length(long.x)
        gic <- log(deviance(fit)/n) + (fit$df *c * log(ncol(long.v)))/n #this is still too aggressive
        #pick the best one
        #I think we want zou + the cleanup thing ?

        #Trying the stabilizing selection- this might be a good choice for the last iterations, once the model has been selected.
        subs <- list()
        for(i in 1:50)
        {
          is <- sample(1:length(long.x), floor((0.9 * length(long.x))));
          subs[[i]] <- glmnet::glmnet(x = long.v[is,], y = long.x[is], family = "gaussian", alpha = 1,
                                      intercept = FALSE, standardize = TRUE,nlambda = 100, trace.it=1)
        }
        z.counts <- lapply(subs, function(x) x$beta !=0)
        per.setting <- do.call("sum", lapply(subs, function(x) x$beta !=0)) #assuming same lambda each time, dubious
        prob.entries <- Reduce("+", z.counts)/50
        highest.prob.by.col <- apply(prob.entries,1, function(x) max(x))
        pi_thresh = 0.9
        drop.entries <- which(!(highest.prob.by.col > pi_thresh))
        coefs.out <- coef(fit, s = alpha.list[min.index])
        coefs.out[drop.entries] <- 0
        u.ret = matrix(coefs.out, nrow = nrow(X), ncol = ncol(V),byrow = TRUE)
      }

#look at the probability for a given
        min.index <- which.min(bic.list)
        #message("Selected BIC index: ", min.index)
        #print(alpha.list)
        #message(alpha.list[min.index][-1])
        #print(option$bic.var)
        #TODO: check this
        u.ret = matrix(coef(fit, s = alpha.list[min.index])[-1], nrow = nrow(X), ncol = ncol(V),byrow = TRUE)
        #return teh best lambda associated with it
        #give the range for the next run, centered around this find
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

        #ProposeNewSparsityParams(bic.list, fit$lambda, (fit$lambda), 2, n.points = 20) #need to modify this to allow for conditions.....
        #Check: is in the realm of those picked by CV?
        return(list("U" = pm(u.ret),"alpha.sel"=fit$lambda[min.index],
                    "bic"= bic.list[min.index], "sparsity_space"=max(fit$lambda), "total.log.lik" = NA,
                    "penalty" = penalty, "s"=s, "next.upper" = upper.next, "next.lower" = lower.next, "SSE"=sse))


      }
      else{ #just fitting a single instance:
        message("Fitting model")
        pryr::mem_change(fit <- glmnet::glmnet(x = long.v, y = long.x, family = "gaussian", alpha = 1,
                             intercept = FALSE, standardize = option$std_coef,lambda = option$alpha1,trace.it=1))
        message("Fitting complete")
        penalty <- fit$penalty
        u.ret = matrix(coef(fit, s = option$alpha1)[-1], nrow = nrow(X), ncol = ncol(V),byrow = TRUE)
        return(list("U" = u.ret, "sparsity_space"=u.ret, "total.log.lik" = NA, "penalty" = penalty, "s"=s, "SSE"=deviance(fit)))
      }
    } else if( option$regression_method == "penalized")
    {
      if(is.null(formerU))
      {

        #Change penalty factor for V
        fit = penalized::penalized(response = long.x, penalized = as.matrix(long.v), lambda1 = option[['alpha1']],lambda2=0, unpenalized = ~0,
                                   positive = FALSE, standardize = option$std_coef, trace = FALSE)
      }else
      {
        fit = penalized::penalized(response = long.x, penalized =  as.matrix(long.v), lambda1 = option[['alpha1']],lambda2=0, unpenalized = ~0,
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

      if(option$actively_calibrating_sparsity) { max.sparsity <- rowiseMaxSparsity(Matrix::Matrix(long.x), (long.v))}
      return(list("U" = pm(l), "sparsity_space"=max.sparsity, "total.log.lik" = ll, "penalty" = penalty, "s"=s, "SSE"=sse))
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
#'#
  FitUWrapper <- function(X,W,W_c,V, option, prevU = NULL,...)
  {
    #if(option[['ncores']] > 1) #This is not working at all. Can't tell you why. But its not. Need to spend some time debugging at some point.
    #{
    #  U = fitUParallel(X, W,W_c, V, option, formerU = prevU); #preL is by default Null, unless yo specify!
#
    #}else
    #{
      #U = fit_U(X, W, W_c, V, option,formerU = prevU, ...);
      U = FitUGlobal(X, W,W_c, V, option,formerU = prevU, ...);

    #}
    return(U)
  }



#Experimental functions
##Examining possible convergence settings
trackAlternativeConvergences <- function(X,W,W_c, U,V, index, record.data, u.alphamax, v.lambdamax, option_old)
{
  message(dim(U))
  message(dim(V))
  track.modes <- list()
  #simple sum mode
  EPSILON = 1e-4
  track.modes[["simple.sum"]] = record.data$alpha.s[[index]] + record.data$lambda.s[[index]]
  track.modes[["U.norm"]] = norm(U, "F")
  track.modes[["V.norm"]] = norm(V, "F")
  track.modes[["U"]] = U #For memory?
  track.modes[["V"]] = V #For memory?
  #scaled sum mode
  track.modes[["lambda.max"]] <-  v.lambdamax
  track.modes[["alpha.max"]] <- u.alphamax
  track.modes[["alpha.scaled"]] <- record.data$alpha.s[[index]]/u.alphamax
  track.modes[["lambda.scaled"]] <- record.data$lambda.s[[index]]/v.lambdamax
  track.modes[["scaled.sum"]] = record.data$alpha.s[[index]]/u.alphamax + record.data$lambda.s[[index]]/v.lambdamax
            track.modes[["bic.alpha"]] = record.data$bic.a[[index]]
              track.modes[["bic.lambda"]] = record.data$bic.l[[index]]

             option <- option_old
              option$fixed_ubiq <- TRUE
              option$alpha1 <- record.data$alpha.s[[index]]
              option$lambda1 <-  record.data$lambda.s[[index]]
              track.modes[["decomp.obj"]] <- compute_obj(X, W, W_c, U, V, option, decomp = TRUE)
              track.modes
}


UpdateAndCheckSparsityParam <-  function(prev.list, new, errorReport = FALSE)
{
 if(!is.null(prev.list))
    {

	 prev <- prev.list[length(prev.list)]
	  new <- new[1]
	if( errorReport & !is.null(prev.list) > 1 & (abs(new-prev)/prev > 1.5))
	  {
	    message("WARNING: large jumpy in sparsity parameter encountered....")
	    message("Current: ", prev)
	    message("New: ", new)
	  }
	}
  return(c(prev.list, new))
}


#Estimator of the L1 degrees of freedom
MatrixDFV <- function(mat_in, fixed_first = FALSE){
  print("RED ALERT- USING THE BAD V")
  if(fixed_first)  #if(FALSE) FF$
  {
    sum(mat_in[,-1] == 0)
  }else
  {
    sum(mat_in == 0)
  }

}
MatrixDFU <- function(mat_in,fixed_first = FALSE)
{
  #this means we are learning U
  #Use the degrees of freedom: number of non-zero coefficients
  #sum(round(mat_in, digits = 4) != 0)
  #if(fixed_first) FF$
  if(fixed_first)  #if(FALSE) FF$
  {
	  sum(mat_in[,-1] != 0)
  }
  else
  {
	  sum(mat_in != 0)
  }
  #This seems like it would favor sparsity?
}

#CalcMatrixBIC(X,W,U,V,df=df,fixed_first = fixed_first,...)
#Note- we assume input to already be in the correct direction, so no transfofmration eeded on the covariance term...
CalcMatrixBIC <- function(X,W,U,V, W_cov = NULL, ev="std", weighted = FALSE, df = NULL, fixed_first = FALSE, learning = "V")
{
  `%>%` <- magrittr::`%>%`
  #Calculated based on my understanding of what is given by equations 12 and 7 in Lee et al 2010
  #We assume U (First term) is known/fixed, V is our predictor we evaluate here
  #11/17 UPDATE: df is calculated for U and V!
  if(!weighted)
  {
    message("unweighted, alert.")
  }
  n=nrow(X)
  d = ncol(X)
  #if(fixed_first) FF$
  #This shouldn't make a difference, resiudals will be unchanged.
  #TODO: remove this.
  if(fixed_first)  #if(FALSE) FF$
  {
	  #message("Removing first factor because its fixed")
	  #Regress out of X the first factor effects.
	  #Remove the first factor from the downstream steps
	  X <- X -  (U[,1] %*% t(V[,1]))
	  V <- V[,-1]
	  U <- U[,-1]
  }

  #its possible the matrices are different sizes, depending on which stage of the fitting they get passed in. Correct for this:
  #aligned <- AlignFactorMatrices(X,W,U,V); U <- aligned$U; V<- aligned$V
  k = ncol(U)
  if(!weighted)
  {
    message("not weighting here..")
    W <- matrix(1, nrow = nrow(X), ncol = ncol(X))
  }
  if(ev =="std" )
  {
    errvar <- AltVar(X,U)
  }else if(ev == "ols" & weighted)
  {
    errvar <- WeightedAltVar(X,W,U, var.method = bic.var, W_cov = W_cov)
  }
  else # anuumber gets passed in.
  {
    message("in correct place...")
    #print(ev)
    errvar = ev
  }

  if(is.null(df))
  {
	  message("something went wrong,..")
    df <- MatrixDFU(V)
  }

  if(is.null(W_cov))
  {
    W_cov <- diag(nrow(X))
  }

  #I was doing this wrong previously.
  #If we reotated things, we need to rotate them back....
  if(learning == "U")
  {
    u.old <- U;
    U <- V ; V <- u.old
    X <- t(X); W <- t(W);
  }
  if(learning == "V")
  {
    message("Currently no covariance on U")
    W_cov = NULL
  }
  resid.matrix <- calcGlobalResiduals(X,W,U,V, W_cov = W_cov, fixed_first = fixed_first)
  stopifnot(any(!is.null(resid.matrix)))
  ret <- norm(resid.matrix, type = "F")^2/(n*d * errvar) + (log(n*d)/(n*d))*df
  #message("Total ", ret)
 #message(". ")

  return(ret)
}

CalcMatrixBIC.loglikversion <- function(X,W,U,V, which.learning, W_cov = NULL, df = NULL, lm.fit.residuals = NULL,...)
{
  if(is.null(W_cov) & which.learning == "U")
  {
    message("Error case that should never happen -giving U but no covar matrix. at minimum identity matrix.")
  }
  `%>%` <- magrittr::`%>%`
  n=nrow(X)
  d = ncol(X)
  k = ncol(U)
  if(which.learning == "V")
  {
    model.ll = penalizedLogLikV(X,W,U,V,...) #pass in the fixed first. Only effect is on the residuals.
  }else
  {
    #learning U
    model.ll = penalizedLogLikU(X,W,W_cov,U,V)
  }
  if(is.null(df))
  {
    df <- MatrixDFU(V)
  }
  #ret <- -model.ll + (log(n*d)/(n*d))*df
  ret <- -2*model.ll + (log(n*d))*df
  if(!is.null(lm.fit.residuals))
  {
    if(which.learning == "V")
    {
      lm.ll = penalizedLogLikV(X,W,U,V, use.resid = lm.fit.residuals,...)
    }else
    {
      #learning U
      lm.ll = penalizedLogLikU(X,W,W_cov,U,V, use.resid = lm.fit.residuals)
    }
    #ret <- neg.ll/(lm.ll) + (log(n*d)/(n*d))*df
    #the goal is to get a ratio which is SMALLER if the fit is better, larger if the fit is worse.
    #Since bigger is better with ll, we actually have to invert it.
    ret <- lm.ll/(model.ll) + (log(n*d)/(n*d))*df
    message("Warning- this BIC method doesn't work well if there are negatives involved...")
    if(model.ll < 0)
    {
        message('havent figured out this case yet, scores will be uninterpretable')
    }
  }
  return(ret)
}


#' Wrapper for all BIC methods under consideration.
#'
#' @param fit glmnet object
#' @param covars long.v or long.u
#' @param y long.x- the explanatory data (weighted)
#' @param bic.mode which version to use; options are Zou, dev, sklearn, or the extended variants of them.
#' @return BIC score for each call (list)
#' @export
calculateBIC <- function(fit, covars,y, bic.mode)
{
    switch(
    bic.mode,
    "sklearn_eBIC"= sklearnBIC(fit,covars,y) + extendedBIC(fit),
    "sklearn"= sklearnBIC(fit,covars,y),
    "Zou_eBIC"= ZouBIC(fit, covars, y) + extendedBIC(fit, scale = TRUE),
    "Zou"= ZouBIC(fit, covars, y),
    "dev_eBIC"= BICglm(fit)+ extendedBIC(fit),
    "dev"= BICglm(fit),
    "std"= stdBIC(fit,covars,y),
    sklearnBIC(fit,covars,y)
          )
}

#' Calculate the extended term which accounts for model size, based on work by Chen et al (2008)
#'
#' @param fit glmnet object
#'
#' @return extended
#' @export
extendedBIC <- function(fit, scale = FALSE)
{
    #p is the number of parameters under consideration
    #n is the
    p = dim(fit$beta)[1] #all the possible parameters to learn
    n = nobs(fit)
    gamma = 1-1/(2 * log(p) / log(n))
    r <- 2 * gamma * lchoose(p, fit$df)
    ifelse(scale, r/nobs(fit), r)
}

#' Calculate the BIC using a vanilla, textbook standard formulation from a GLMNET fit object
#'
#' @param fit the GLMNET fit object
#' @param long.v the covariates matrix, sparse
#' @param long.x the outcome variable
#'
#' @return a vector of BIC scores across multiple different \lambda settings
#' @export
#'
stdBIC <- function(fit, long.v, long.x)
{
  all.preds <- predict(fit, newx = long.v)
  #I've already extended the dims here, so n = MN
  p = 1
  n = nobs(fit)
  sapply(1:ncol(all.preds), function(i) -2*penLLSimp(n, p, long.x - all.preds[,i]) + log(n) * fit$df[i])
}

#' BIc calculation for use with glmnet
#'
#' @param fit glmnet object
#'
#' @return BIC score for each call
#' @export
BICglm <- function(fit, bic.mode = ""){
  #based on code from https://stackoverflow.com/questions/40920051/r-getting-aic-bic-likelihood-from-glmnet
  tLL <- -deviance(fit) # 2*log-likelihood
  k <- fit$df
  n <- nobs(fit)
  BIC <- log(n)*k - tLL
  BIC
}

# sklearnBIC(fit, long.v, long.x)

#' Calculate BIC based on what is in the sci-kit learn package (see https://scikit-learn.org/stable/modules/linear_model.html#lasso-lars-ic)
#' this differs from the traditional derivation, in that the use the unbiased estimate for s rather than the MLE.
#' Source code is here https://github.com/scikit-learn/scikit-learn/blob/872124551/sklearn/linear_model/_least_angle.py#L1991
#'
#' @param fit
#' @param x
#' @param y
#' @param bic.mode
#'
#' @return
#' @export
#'
#' @examples
sklearnBIC <- function(fit,x,y, bic.mode = "")
{
  p = ncol(x)
  n = length(y)
  #OLS fit
  #This is expensive to do afeter taking so long to fit a goshdarn lasso. Need to find an alternative option?
  lm.fit <- glmnet::bigGlm(y=y , x=x, family = "gaussian",intercept = FALSE,trace.it = 1)
  fit.best.lm <- predict(lm.fit, newx = x)

  #our fit
  coef = coef(fit)
  yhat=x%*%coef[-1,]
  residuals = (y- yhat)
  sse = apply(residuals, 2, function(x) sum(x^2))


  self.noise.variance <- sum((y - fit.best.lm)^2)/(n -p)
  scikitlearn.bic <- n*log(2*pi*self.noise.variance) + sse/self.noise.variance + log(n)*fit$df
  scikitlearn.bic

}

ZouBIC <- function(fit, x, y, var.meth = "mle")
{
  n <- nobs(fit)
  p=ncol(x)
  #lm.fit <- predict(fit, newx = x, s= 0)
  lm.fit <- glmnet::bigGlm(x, y, family = "gaussian", intercept = FALSE)
  all.preds <- predict(fit, newx = x) #pretty slow step, surprisingly.
  bic.new <- c()
  fit.new <-  c()
  df.term <- c()
  #best.fit <- norm(y-lm.fit, "2")^2/n
  best.fit <- deviance(lm.fit)/n
  if(var.meth == "unbiased")
  {
    best.fit <- deviance(lm.fit)/(n - length(x))
  }
  for(i in 1:length(fit$df))
  {
    curr.fit <- norm(y-all.preds[,i], type ="2")^2/n
    left.term <- curr.fit / best.fit
    fit.new <- c(fit.new, left.term)
    right.term <- (log(n) * fit$df[i])/n
    df.term <- c(df.term, right.term)
    bic.new = c(bic.new, left.term + right.term)
  }
  bic.new
}


#' BIC using global fit calculation instead of sum of each regression's liklihood
#' Note that this doesn't require any fancy counting on the residuals, since all are just lumped into one regression
#' @param X original data matrix (SNPs x studies)
#' @param W uncertainty associated with X (SNPs x studies)
#' @param U SNPs x K matrix
#' @param V Studies x K matrix
#' @param W_cov Covariance matrix to correct for. For now just implemented in U.
#' @param which.learning if we are learning "V" or "U"
#' @param df degrees of freedom associated with matrix currently learning
#' @param lm.fit.residuals Specify these to scale the log-liklihood by the optimal OLS fit. Not recommended.
#'
#' @return
#' @export
CalcMatrixBIC.loglikGLOBALversion <- function(X,W,U,V, which.learning, W_cov = NULL, df = NULL, lm.fit.residuals = NULL, decomp = FALSE,...)
{
  `%>%` <- magrittr::`%>%`
  n=nrow(X)
  d = ncol(X)
  k = ncol(U)
  #Residuals are the same regardless of
  #extend_the_weighting
  #remove conditions below
  resids = calcGlobalResiduals(X,W,U,V,W_cov=W_cov,...)

  #if(which.learning == "U")
  #{
  #  resids = calcGlobalResiduals(X,W,U,V,W_cov=W_cov,...)
  #}else
  #{
    #why is this condition here? might just make thigns complex.
  #  resids = calcGlobalResiduals(X,W,U,V, W_cov = NULL,...)
  #}

  #model.ll = penLL(n*d, resids)
  model.ll = penLLSimp(n,d,resids)
  #The following was tested on 3/17. It doesn't work, my derivation is wrong. I am missing something in the log.
  #message("Compare this to...")
  #print(penLLEmpirical(n*d, resids))
  #message("Trying to look at things differently")
  #model.ll=penLLEmpirical(n*d, resids)
  if(is.null(df))
  {
    df <- MatrixDFU(V)
  }
  #message("Not scaling by n*d here anymore...")
  #ret <- -model.ll + (log(n*d)/(n*d))*df
  if(decomp)
  {
    message("log fit: ", -2*model.ll)
    message("df term: ", log(n*d)*df)
  }
  if(which.learning == "U" & df == 0)
  {
    message("Red alert: we have zeroed out everything")
    message("This means our model thinks there is no signal whatsoever. Also possible.")
  }
    ret <- -2*model.ll + log(n*d)*df
    if(!is.null(lm.fit.residuals))
    {
      lm.ll = penLL(n*d, lm.fit.residuals)
      ret <- lm.ll/(model.ll) + (log(n*d)/(n*d))*df
    }
  return(ret)
}



#' BIC using global fit calculation and the optimizing log-liklihood (a normal form)
#' Note that this doesn't require any fancy counting on the residuals, since all are just lumped into one regression
#' @param X original data matrix (SNPs x studies)
#' @param W uncertainty associated with X (SNPs x studies)
#' @param U SNPs x K matrix
#' @param V Studies x K matrix
#' @param W_cov Covariance matrix to correct for. For now just implemented in U.
#' @param which.learning if we are learning "V" or "U"
#' @param df degrees of freedom associated with matrix currently learning
#' @param lm.fit.residuals Specify these to scale the log-liklihood by the optimal OLS fit. Not recommended.
#'
#' @return
#' @export
CalcMatrixBIC.NormalLoglikGLOBALversion <- function(X,W,U,V, which.learning, W_cov = NULL, df = NULL, calc.residual.variance = FALSE, decomp = FALSE,...)
{
  `%>%` <- magrittr::`%>%`
  n=nrow(X)
  d = ncol(X)
  k = ncol(U)
  #Residuals are the same regardless of
  if(which.learning == "U")
  {
    resids = calcGlobalResiduals(X,W,U,V,W_cov=W_cov,...)
  }else
  {
    resids = calcGlobalResiduals(X,W,U,V, W_cov = NULL,...)
  }
  residual.variance =1
  if(calc.residual.variance)
  {
    residual.variance = CalcVariance(n*d, resids)
  }

  model.ll = stdLogNormalFit(resids, residual.variance=residual.variance)
  if(is.null(df))
  {
    df <- MatrixDFU(V)
  }

  if(decomp)
  {
    message("log fit: ", -2*model.ll)
    message("df term: ", log(n*d)*df)
  }
  if(which.learning == "U" & df == 0)
  {
    message("Red alert: we have zeroed out everything")
    message("This means our model thinks there is no signal whatsoever. Also possible.")
  }
  ret <- -2*model.ll + log(n*d)*df
  if(!is.null(lm.fit.residuals))
  {
    lm.ll = penLL(n*d, lm.fit.residuals)
    ret <- lm.ll/(model.ll) + (log(n*d)/(n*d))*df
  }
  return(ret)
}




#Determine the variance using OLS. Basically, this is the "best" residual variance we can get.
AltVar <- function(X,U)
{
  resids <- c()
  for(col in 1:ncol(X))
  {
    fit <- lm(X[,col]~U)
    resids <- c(resids, resid(fit))
  }
  n = length(resids)
  (var(resids)* (n - 1)) / n
}

#WeightedAltVar(t(X),t(W),initV)
#Determine the variance using OLS. Basically, this is the "best" residual variance we can get.

WeightColumnHelper <- function(w,x,U)
{
  return(list("wx"=w*x, "wu" = w*U))
}

#Meant to help regress out one learned component at a time.
#not actually that helpful
#mc is the matrix component to pull out.
RegressOutMatrixComponent <- function(X,W,mc)
{
  ret <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
  for(col in 1:ncol(X))
  {
    w <- unlist(W[,col])
    wx <- (w*X[,col])
    wu <- w * mc
    fit <- lm(wx~wu + 0)
    fitted <- mc*coef(fit)
    ret[,col] <- X[,col] - fitted
  }
  return(ret)
}

WeightedAltVar <- function(X,W,U, var.method = "mle", fit.resids = FALSE, W_cov = NULL)
{
  if(is.null(W_cov))
  {
    W_cov <- diag(nrow(X))
    message("filling in with blank W")
    print(dim(W_cov))
  }
  n <- nrow(X) * ncol(X)
  resids <- c()
  fit.mat.dat <- NULL
  for(col in 1:ncol(X))
  {
    w <- unlist(W[,col])
    wx <- W_cov %*% (w*X[,col])
    wu <- W_cov %*% (w * U)
    fit <- lm(wx~wu + 0)
    resids <- c(resids, resid(fit))
    if(fit.resids)
    {
      fit.mat.dat <- rbind(fit.mat.dat, unlist(residuals(fit)))
    }
  }
  stopifnot(length(resids) == n)
  r = 0
  if(var.method == "unbiased")
  {
    message("unbiased")
    p <- ncol(U) * ncol(X)
    r = (1/(n-p))*sum(resids * resids)
  }else if(var.method == "map")
  {
    message("map")
    alpha <- 0.001
    beta <- 0.001
    den <- (alpha + n/2 - 1)
    num <- beta + 0.5*sum(resids * resids)
    r <- num/den
  }else if(var.method == "mle") #MLE
  {
    #message("mle")
    p <- 0
    r = sum(resids * resids)/n
  } else #var method = std
  {
    message("std method.")
    r <- var(resids)
  }

  #print("my calculated r is:")
  #print(r)
  #print("mle var")
  #print((var(resids) * (n-1) )/ n)
  #print(r)
  #This is equivilant to all ther terms being put into one matrix
  if((sum(resids * resids) == 0))
  {
    message("No residual variance with current fit.")
    print(r)
    message("Be wary here- not clear what to do. Going to just approximate to small number")
    message("This is strong evidence of overfitting- we recommend dropping current factors.")
    message("Current U dims:")
    print(dim(U))
    #this would be easy enough to do, if we had everything in "PVE" order.
    r = 1e-10
  }
  if(fit.resids)
  {
    return(list("resid.var" = r, "resids" = as.matrix(fit.mat.dat)))
  }
  return(r)
}

#Helper calls
CalcVBIC <- function(X,W,U,V,fixed_first=FALSE,lm.resid = NULL,...)
{
  #message("Harding coding fixed_first = FALSE on all BIC calculations")
  #fixed_first = FALSE
  df.dat=MatrixDFU(V,fixed_first=fixed_first)
  #var.based.bic <- CalcMatrixBIC(X,W,U,V,df=df.dat,fixed_first = fixed_first,learning = "V",...)
  var.based.bic <- NA
  #ll.based.bic <- CalcMatrixBIC.loglikversion(X,W,U,V, which.learning = "V", df = df.dat, fixed_first=fixed_first)

  ll.based.bic <- NA
  #ll.based.bic.ratio <- CalcMatrixBIC.loglikversion(X,W,U,V, which.learning = "V", df = df.dat, lm.fit.residuals = lm.resid,fixed_first=fixed_first)
  ll.based.bic.ratio <- NA
  ll.global.bic <- CalcMatrixBIC.loglikGLOBALversion(X,W,U,V, which.learning = "V", df = df.dat,fixed_first=fixed_first)
  #ll.global.ratio.bic <- CalcMatrixBIC.loglikGLOBALversion(X,W,U,V, which.learning = "V", df = df.dat,
                                                           #lm.fit.residuals = lm.resid,fixed_first=fixed_first) #OMIT THIS.
ll.global.ratio.bic <- NA
  return(list(var.based.bic,ll.based.bic,ll.based.bic.ratio,ll.global.bic,ll.global.ratio.bic)) #WE LIKE options 1,2,4
}

CalcUBIC <- function(X,W,W_c, U,V,lm.resid = NULL,...)
{
  df.dat <-  MatrixDFU(U)
  #var.based.bic <- CalcMatrixBIC(t(X),t(W),V,U,df=df.dat,W_cov = W_c, learning = "U",...)
  #And now, the alternative versions:
  #ll.based.bic <- CalcMatrixBIC.loglikversion(X,W,U,V, which.learning = "U", df = df.dat, W_cov = W_c)
  #ll.based.bic.ratio <- CalcMatrixBIC.loglikversion(X,W,U,V, which.learning = "U", df = df.dat, lm.fit.residuals = lm.resid, W_cov = W_c)#OMIT THIS
  ll.global.bic <- CalcMatrixBIC.loglikGLOBALversion(X,W,U,V, which.learning = "U", df = df.dat, W_cov = W_c)
  #ll.global.ratio.bic <- CalcMatrixBIC.loglikGLOBALversion(X,W,U,V, which.learning = "U", df = df.dat,lm.fit.residuals = lm.resid, W_cov = W_c) #OMIT THIS.

  var.based.bic <- NA
  ll.based.bic <- NA
  ll.based.bic.ratio <- NA
  ll.global.ratio.bic <- NA
  return(list(var.based.bic,ll.based.bic,ll.based.bic.ratio,ll.global.bic,ll.global.ratio.bic)) #WE LIKE options 1,2,4
}


#Fit the V with the scheme:
#Iniital estimates from burn.in.sparsity and consider.parsm
FitVs <- function(X, W,W_c, initU, lambdas,option, weighted = FALSE,reg.elements=NULL)
{
  bic.var = option$bic.var
  f.fits <- list()
  bics <- c()

  for(i in 1:length(lambdas))
  {
    l <- lambdas[[i]]
    option$lambda1 <- l
    # f.fits[[i]] <- fit_V(X, W, initU, option, formerV = NULL) #I forgot to change this: need to re-run objective tests now >_<
    f.fits[[i]] <- FitVWrapper(X, W,W_c, initU, option, formerV = NULL,reg.elements=reg.elements)
  }

  if(weighted) {
    #option$fixed_ubiq <- FALSE
    #if(option$fixed_ubiq) FF$
    message("CURRENTLY SKIPPING the fixed ubiq accomodations for BIC calculation.")
    if(option$fixed_ubiqs & FALSE)  #if(FALSE) FF$
    {
      #If we are down to just 1 column, don't calculate BICs, there is no point anymore.
      if(ncol(initU) == 1)
      {
        #If we have fixed_ubiq and down to 1 column, it doesn't matter
        message("Down to just 1 column, all BIC the same now.")
        return(list("fits" = f.fits, "BIC"=rep(0,length(lambdas))))
      }
      Xr <- RegressOutMatrixComponent(X,W,initU[,1])
      #this should be different for each one, because the fixed column differs
      #this would unfairly advantage fits with more information in factor 1. Doesn't work.
      #In practice, I think this is not different at all from just calculating the weighted variance normally, except
      #That in the variance calculation, they no longer get the benefit of that first factor
      #message("Note: these functions may need to be need to be adjusted- each V1 is different, so each Xr is different that would be learned.")
      #Current setup seems reasonable, but best would be to have 1 for each
      av <- WeightedAltVar(Xr,W,as.matrix(initU[,-1]), var.method = bic.var, fit.resids  = TRUE)
    }else
    {
      av <- WeightedAltVar(X,W,initU, var.method = bic.var, fit.resids  = TRUE)
    }
    #bics <- do.call("rbind", lapply(f.fits, function(x) CalcVBIC(X,W,initU,x$V, ev=av[[1]], lm.resid = av[[2]], weighted = TRUE, fixed_first = option$fixed_ubiqs)))
    bics <- do.call("rbind", lapply(f.fits, function(x) CalcVBIC(X,W,initU,x$V/x$s, ev=1, lm.resid = av[[2]], weighted = TRUE, fixed_first = option$fixed_ubiqs)))
  }else{ #unweighted
    av <- AltVar(X,initU)
    bics <- do.call("rbind", lapply(f.fits, function(x) CalcVBIC(X,W,initU, x$V/x$s, ev=av[[1]],lm.resid = av[[2]], fixed_first = option$fixed_ubiqs)))
  }
  if(Inf %in% bics)
  {
    message("Error here: INF in BIC")
    quit()
  }
  return(list("fits" = f.fits, "BIC"=bics))
}
#Helper function to get rid of empty columns in the data.
#Empty means all the terms are 0.
DropEmptyColumns <- function(matin)
{
  matin <- as.matrix(matin)
  c <- colSums(matin != 0)
  if(any(c == 0))
  {
    #print("dropping")
    drops <- which(c == 0)
    if(1 %in% drops){
      message("BEWARE- 1st column dropping???")
    }
    return(as.matrix(matin[,-drops]))
  }
  return(matin)
}

#Same as the above, but for magrittr piping (not actually used.)
DropEmptyColumnsPipe <- function(lin)
{
  ret.lin <- lin
  ret.lin[[1]] <- DropEmptyColumns(lin[[1]])
  return(ret.lin)

}

#Recalculate the sparsity params for U
#U returned is res-caled up
FitUs <- function(X, W, W_c, initV, alphas,option, weighted = FALSE, reg.elements = NULL)
{
  bic.var = option$bic.var
  l.fits <- list()
  for(i in 1:length(alphas))
  {
    a <- alphas[[i]]
    message(i)
    option$alpha1 <- a
    #l.fits[[i]] <- fit_U(X, W, W_c, initV, option)
    #change made here...3/8
    l.fits[[i]] <- FitUWrapper(X, W, W_c, initV, option,reg.elements=reg.elements)

  }
  #TODO: recode this, so don't need the logic statement. Downstream should be able to handle it
  if(weighted) {
    av <- WeightedAltVar(t(X),t(W),initV, var.method = bic.var, fit.resids  = TRUE, W_cov = W_c)
    #bics <- do.call('rbind', lapply(l.fits, function(x) CalcUBIC(X,W,W_c,as.matrix(x$U),initV, ev=av[[1]], weighted = TRUE, lm.resid=av[[2]])))
    bics <- do.call('rbind', lapply(l.fits, function(x) CalcUBIC(X,W,W_c,as.matrix(x$U / x$s),initV, ev=1, weighted = TRUE, lm.resid=av[[2]])))
  }else{
    av <- AltVar(t(X),initV, fit.resids = TRUE) #This step is quite slow.... need to speed this up somehow.
    bics <-  do.call('rbind', lapply(l.fits, function(x) CalcUBIC(X,W,W_c,as.matrix(x$U/ x$s),initV,ev=av[[1]], lm.resid = av[[2]] )))
  }

  return(list("fits" = l.fits, "BIC"=bics, "resid_var" =av))
}

#From new distribution and current list, how to pick the new ones?
#Bic. list: list of BIc scores for all choices
#optimal.sparsity.param- the top parameter chosen
#new.dist: distribution of all the ne lambda parameter space
#@param curr.mode= the mode of the sparsity space based on the current V and U settings
#@return a list of new sparsity points to try out.


ProposeNewSparsityParams <- function(bic.list,sparsity.params, curr.dist, curr.iter, n.points = 7, no.score = FALSE, one.SD.rule = FALSE, no.drop = FALSE)
{
  curr.mode = DensityMode(curr.dist)
  global.min <- min(curr.dist)
  if(length(bic.list) == 1)
  {
    message("No list to choose from- have zeroed all out..?")
    #Go from cuyrrent value to the mode, give a spread
    return(sort(10^seq(log10(sparsity.params),log10(curr.mode),length.out=n.points)))
  }
  if(no.score)
  {
    #then bic.list is the optimal one; generate fake scores
    fake.scores <- rep(100,length(sparsity.params))
    fake.scores[which(sparsity.params == bic.list)] <- -1
    bic.list <- fake.scores
  }
  optimal.index <- selectOptimalScoreIndex(bic.list, sparsity.params, one.SD.rule)
  #cases with redundancy are complex.
  optimal.sparsity.param <- sparsity.params[optimal.index]
  sorted.sparsity.params <- sort(sparsity.params, index.return = TRUE)
  ordered.list <- sorted.sparsity.params$x
  sorted.index <- which(sorted.sparsity.params$ix == optimal.index)
  #what is the index int eh sorted list of my optimal sparsty parameter?
  #New paradigm: always look above and below,
  #If its the smallest paramter tested
  if(min(ordered.list) == optimal.sparsity.param)
  {
    message('best case is the minimum..')
    if(no.drop & curr.iter > 2)
    {
      message("Not allowing a decrease in score if we are after the first iteration")
      message("Method needs to deal with the minimum as the option")
      above <- ordered.list[sorted.index + 1]
      below <- optimal.sparsity.param
    } else
    {
      above <- ordered.list[sorted.index + 1]
      #below <- optimal.sparsity.param - (above - optimal.sparsity.param)
      #simplify this: we are stil searching, so look orders of magnitude
      #below <- 1e-10
      below <- global.min #maybe a better way to do this
      if(below > above)
      {
        message("Global min param of distribution is greater than current one.")
        message("Setting new minimum to 0.1 of current parameter")
        #print(above)
        #print(below)
        below <- optimal.sparsity.param * 0.1
      }
    }

  } else if(max(ordered.list) == optimal.sparsity.param) #its the largest parameter tested
  {
    message('best case is the maximum')
    below <- ordered.list[sorted.index - 1]
    above <- curr.mode
    #new.list <- 10^seq(log10(below),log10(above),length.out=n.points)
  }else {
    #Its bounded- our estimates should be between the one immediately above and below
    above <- ordered.list[sorted.index + 1]
    below <- ordered.list[sorted.index - 1]
    #new.list <- seq(below,above,length.out = n.points)
  }
  if(length(above) > 1 | length(below) < 1)
  {
    message("WHAT is going on...")
    print(above)
    print(below)
    print(bic.list)
    print(sparsity.params)
    readline()
    quit()
  }
  if(is.na(above) | is.na(below))
  {
    print("proposed new paramters are not possible")
    quit()
  }
  if(above == below)
  {
    message("Converged on single solution")
    return(below)
  }
  new.list <- 10^seq(log10(below),log10(above),length.out=n.points)

  #Ensure none of them are less than 0


  if(any(new.list <= 0))
  {
    rep.list <- new.list[which(new.list > 0)]
    if(any(new.list == 0))
    {
      rep.list <- c(min(rep.list)/2,rep.list)
    }
    new.list <- rep.list
  }

  unique(sort(c(optimal.sparsity.param, new.list)))
}


    oneSDRule <- function(bics, params)
    {
      if(Inf %in% params)
      {
        message("issue here....")
      }
      if(length(unique(bics)) == 1)
      {
        #all the parameters are the same
        message("BIC score for all parameters are same. Likely that all factors have been 0'd out")
        message("Seek recourse.")
        return(which.min(params))
      }
        sd <- sd(bics)
        opt <- min(bics)
        in.range <-bics[(bics > (opt - sd)) & (bics < (opt+sd))]
        if(FALSE)
        {
          print("Optimal")
          print(opt)
          print("SD")
          print(sd)
          print(in.range)
          #top.indices <- which(bics %in% in.range)
          print("all BICs")
          print(bics)
          print("tops")
          print(top.indices)
          print(params)
          #optimal.l <- max(params[top.indices])
          print("Selecting params:")
          print(optimal.l)
        }
        top.indices <- which(bics %in% in.range)
        optimal.l <- max(params[top.indices])
        return(which(params == optimal.l))
    }
#Deals with cases if redundant scores.
#If these are at the upper extreme of the parameter list (likely occurs when all terms have been zeroed out), pick the SMALLEST parameter
#if these are at the lower extreme of the parameter list (likely occurs when the terms are fully dense), pick the largest parameter
  SelectBICFromIdenticalScores <- function(bic, params)
    {
      best.score <- min(bic)
      optimal.index <- which.min(bic)
      #message("Warning- BIC scores are unchanging for certain settings. This is likely evidence of no sparsity.")

      #Choose next lowest bic score index
      i <- sort(bic, index.return = TRUE)

      sorted.bic <- bic[i$ix]
      params.sorted.by.bic <- params[i$ix]
      matching.score.indices <- which(sorted.bic == best.score)

      #If the scores are on the  bigger end of scale
      #if all scores yield the same, its not obvious if we have 0d out or total density. Pick the one closest ot he averagge
      if(length(unique(bic)) == 1)
      {
        message("All scores yield the same BIC. Unclear what to do...")
        mid = abs(params - mean(params))
        optimal.param <- params[which(min(mid) == mid)][1]
      }
      else if(all(min(params.sorted.by.bic[matching.score.indices]) > params.sorted.by.bic[-matching.score.indices]))
      {
        message("Suspect that scores are zeroing out the results, picking the smallest parameter with low BIC")
        optimal.param <- min(params.sorted.by.bic[matching.score.indices])
      } else if(all(max(params.sorted.by.bic[matching.score.indices]) < params.sorted.by.bic[-matching.score.indices]))
      {
        message("Suspect that scores are inducing no sparsity, picking the largest parameter with low BIC")
        optimal.param <- max(params.sorted.by.bic[matching.score.indices])
      }
      else
      {
        #weird case, in the middle. #This means that all of them are equally goood?
        #I this case, I want to pick the most spare one actually
       #print("Beware, unusual case...")
        #optimal.param <- params.sorted.by.bic[matching.score.indices][ceiling(length(matching.score.indices)/2)] #get the middle one
        optimal.param <- max(params.sorted.by.bic[matching.score.indices])
        #print(sorted.bic)
        #print(optimal.param)
        #print(params.sorted.by.bic)
        #print("")
        message("Unusual case with the center, pick the middle. Likely swung too far above or below. Hope is lost :(")

      }

  which(params == optimal.param) #return the optimal index

    }

    #This gets the index for the optimal
  selectOptimalScoreIndex <- function(bic, params, oneSD, ndigits = 6)
  {
    bic <- round(bic, digits = ndigits)
    if(oneSD)#just testing our the one sd rule
    {
      optimal.index <- oneSDRule(bic,params)
    } else
    {
      best.score <- min(bic)
      optimal.index <- which.min(bic)
      if(length(which(bic == best.score)) > 1)
      {
        optimal.index <- SelectBICFromIdenticalScores(bic, params)
        if(length(optimal.index) > 1)
        {
          message("Warning: proposing redundant values.")
          optimal.index <- optimal.index[1]
        }
      }
    }
    optimal.index
  }
    #This function selects the optimal matrix and drops non-zero entries.
  selectOptimalInstance <- function(fit.data, bic, params, oneSD = FALSE)
  {

    optimal.index <- selectOptimalScoreIndex(bic, params, oneSD)
    if(all(fit.data$fits[[optimal.index]]$s) == 1)
    {
      optimal.matrix <- DropEmptyColumns(fit.data$fits[[optimal.index]][[1]])
    }else
    {
      which.to.drop <- which(colSums(fit.data$fits[[optimal.index]][[1]] != 0) == 0)
      optimal.matrix <- DropEmptyColumns(fit.data$fits[[optimal.index]][[1]]) / fit.data$fits[[optimal.index]]$s[,-which.to.drop]
    }
    if(ncol(optimal.matrix) == 0)
    {
      optimal.matrix <- matrix(0, nrow = nrow(fit.data$fits[[optimal.index]][[1]]), ncol = 1 )
    }
    #SELECT NON-ZERO entries
  return(list("m" = optimal.matrix, "p" = params[[optimal.index]], "index" = optimal.index))
  }


#  bic.dat <- getBICMatrices(opath,option,X,W,all_ids, names)
  #getBICMatrices(opath,option,X,W,all_ids, names, burn.in.iter = 1)
  #option$bic.var <- "mle". #unbiased is oto strong here
getBICMatrices <- function(opath,option,X,W,W_c, all_ids, names, min.iter = 2,
                           max.iter = 100, burn.in.iter = 0, bic.type = 4, use.optim = TRUE, rg_ref = NULL,reg.elements = NULL)
{
  message("using BIC type:  ", bic.type)
  W_ld = NULL
  conv.options <- list()
#If we get columns with NA at this stage, we want to reset and drop those columns at the beginning.
  #Currently, just using NULL for W_ld
  burn.in.sparsity <- DefineSparsitySpaceInit(X, W, W_c, NULL, option, burn.in = burn.in.iter, rg_ref =rg_ref,reg.elements=reg.elements ) #If this finds one with NA, cut them out here, and reset K; we want to check pve here too.

  #optimal.v <- DropLowPVE(X,W,burn.in.sparsity$V_burn, thresh = 0.01)  #skipping... in sims seems to be bad?
  if(option$u_init != "")
  {
    optimal.u <- burn.in.sparsity$U_burn
    option$K <- burn.in.sparsity$new.k
    message("initializing with U")
  }else
  {
    optimal.v <- burn.in.sparsity$V_burn
    #Doesn't appear the order makes a big difference.
    #u.sparsity <- DefineSparsitySpace(X,W,as.matrix(optimal.v),"U", option)
    option$K <- ncol(optimal.v)
  }
  consider.params <- SelectCoarseSparsityParams(burn.in.sparsity, burn.in.iter, n.points = 15)

  #things to record
  rec.dat <- list("alphas"=c(), "lambdas"=c(), "bic.a" = c(), "bic.l"=c(), "obj"=c(),
                  "v.sparsity" = c(), "u.sparsity"=c(), "iter"=c(), "sd.sparsity.u" = c(), "sd.sparsity.v" = c(),
                  "alpha.s" = c(), "lambda.s" = c(), "Ks" = c(), "Vs" = list(), "X_hat" = c(), "sparsity.obj" <- list())
  #kick things off
  lambdas <- consider.params$lambdas
  alphas <- consider.params$alphas
  INCR.LIMIT=3
  OPTIMIZER="SANN"

  #DECR.LIMIT=4 Think about ths some more...
  NOT.CONVERGED <- TRUE; i = 1
  #CONV.MODE = "BIC.change"
  CONV.MODE <- option$param_conv_criteria
  message("Convergence mode: ", CONV.MODE)
  #If initializing with U, start there...
  if(option$u_init != "")
  {
    rec.dat$lambdas <- c(rec.dat$lambdas, lambdas); #save current alphas.
    v.fits <- FitVs(X,W,W_c, optimal.u,lambdas,option, weighted = TRUE)
    #Which versin of BIC do we want to use this time?
    optimal.iter.dat <- selectOptimalInstance(v.fits, unlist(v.fits$BIC[,bic.type]), lambdas)
    optimal.v <- optimal.iter.dat$m
    rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v))
    rec.dat$lambda.s <- c(rec.dat$lambda.s,optimal.iter.dat$p)

    #Record new data
    rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
    rec.dat$bic.l <- c(rec.dat$bic.l,unlist(v.fits$BIC[,bic.type]))
    #update the parameters for U based on the new V
    u.sparsity <- DefineSparsitySpace(X,W,W_c, optimal.v, "U", option,reg.elements=reg.elements) #Here in case we hit the max.
    alphas <- SelectCoarseSparsityParamsGlobal(u.sparsity, n.points = 15)
  }
  #$Remove low PVE right now.
  #no.drop.live = TRUE
  no.drop.live = FALSE
  #message("Swapped 15 points for 15 points for speed.")
  while(i < max.iter & NOT.CONVERGED){
    print("iter on now...")
    print(i)
    if(i == 22)
    {
      message('getting down to it....')
    }
    #now fit U:
    rec.dat$Vs[[i]] <- optimal.v

    #Now fit U!
    if(use.optim)
    {
      #u.fits <- FitUs(X, W, optimal.v, alphas,option, weighted = TRUE)
      #bic.list.u <- unlist(u.fits$BIC[,bic.type])
      #^above for debugging and comparison.
      if(i == 1)
      {
        init.alpha = alphas[floor(length(alphas)/2)]
        OPTIMIZER="SANN"
        #OPTIMIZER="Brent"
        #init.alpha = min(alphas)
      }else
      {
        init.alpha = rec.dat$alpha.s[length(rec.dat$alpha.s)] #the most recent previous guess
        OPTIMIZER="Brent"
      }
      upper.lim=max(alphas)
      if(i > 1)
      {
        #Don't allow the parameter to more than triple in a given iteration
        upper.lim = min(init.alpha*INCR.LIMIT, max(alphas))
      }
      #scoretest <- GetUBICOptim(init.alpha, X,W,optimal.v, option )
      #Trying this as SANN instead ob rent- meant to be bette ron rough surfaces.
      #message('here.')
      test <- optim(par = init.alpha, fn =  GetUBICOptim, method = OPTIMIZER, lower = min(alphas)*0.1, upper = upper.lim,
                    X=X, W=W, initV = optimal.v, option = option, W_c = W_c, control= list('trace'=1))
      #test.brent <- optim(par = init.alpha, fn =  GetUBICOptim, method = OPTIMIZER, lower = min(alphas)*0.1, upper = upper.lim,
                    #X=X, W=W, initV = optimal.v, option = option, W_c = W_c, control= list('trace'=1))
      #test.min <- optim(par = min(alphas), fn =  GetUBICOptim, method = OPTIMIZER, lower = min(alphas)*0.1, upper = upper.lim,
                    #X=X, W=W, initV = optimal.v, option = option, W_c = W_c, control= list('trace'=1))
      #test.brent.min <- optim(par = min(alphas), fn =  GetUBICOptim, method = "Brent", lower = min(alphas)*0.1, upper = upper.lim,
                    #X=X, W=W, initV = optimal.v, option = option, W_c = W_c, control= list('trace'=1))
      #test <- optim(par = init.alpha, fn =  GetUBICOptim, method = "Brent", lower = min(alphas)*0.1, upper = upper.lim,
      #              X=X, W=W, initV = optimal.v, option = option, bic.method = 1, ev = 1, W_c = W_c, control= list('trace'=1))
      message("Optimizer is ", OPTIMIZER)
      u.fits <- FitUs(X, W, W_c, optimal.v, c(test$par),option, weighted = TRUE)
      bic.list.u <- c(test$value)
      alphas <- c(test$par)
    }else
    {
      u.fits <- FitUs(X, W, W_c, optimal.v, alphas,option, weighted = TRUE)
      bic.list.u <- unlist(u.fits$BIC[,bic.type])

    }
    rec.dat$alphas <- UpdateAndCheckSparsityParam(rec.dat$alphas, alphas, errorReport = TRUE)

    rec.dat$bic.a <- c(rec.dat$bic.a,bic.list.u)
    #Pick the best choice from here, using threshold.
    optimal.iter.dat <- selectOptimalInstance(u.fits, bic.list.u, alphas)
    optimal.u <- optimal.iter.dat$m
    rec.dat$alpha.s <- c(rec.dat$alpha.s,optimal.iter.dat$p)

    rec.dat$U_sparsities = c(rec.dat$U_sparsities, matrixSparsity(optimal.u, ncol(X)));

    #now get the new lambdas for V:
    #Is this what we want?
    v.sparsity <- DefineSparsitySpace(X,W,W_c, as.matrix(optimal.u),"V", option,reg.elements=reg.elements)
    if((i == 1 & option$u_init == "") | use.optim)
    {
      lambdas <- SelectCoarseSparsityParamsGlobal(v.sparsity, n.points = 15)
    }else
    {
      lambdas <- ProposeNewSparsityParams(bic.list.v, lambdas, v.sparsity, i, n.points = 7, no.drop = no.drop.live)
    }
    rec.dat$sd.sparsity.v <- c(rec.dat$sd.sparsity.v,sd(v.sparsity))

    #Now fit V!
    if(use.optim)
    {
        #v.fits <- FitVs(X,W, optimal.u,lambdas,option, weighted = TRUE)
        #bic.list.v <- unlist(v.fits$BIC[,bic.type])
      if(i == 1)
      {
        init.lambda = lambdas[floor(length(lambdas)/2)]
        OPTIMIZER="SANN"
        #OPTIMIZER="Brent"
        #init.lambda = min(lambdas)
      }else
      {
        OPTIMIZER="Brent"
        init.lambda = rec.dat$lambda.s[length(rec.dat$lambda.s)]
      }
      upper.lim=max(alphas)
      if(i > 1)
      {
        #Don't allow the parameter to more than triple in a given iteration
        upper.lim = min(init.lambda*INCR.LIMIT, max(lambdas))
      }
      message("Optimizer is ", OPTIMIZER)
        #scoretest <- GetVBICOptim(init.lambda, X,W,optimal.u, option )
      #v.fits <- FitVs(X,W, optimal.u,min(lambdas),option, weighted = TRUE)
        #test <- GetVBICOptim(par, X,W,optimal.u, option, weighted = TRUE, bic.method = 4)
        test <- optim(par = init.lambda, fn =  GetVBICOptim, method = OPTIMIZER, lower = min(lambdas)*0.1, upper = upper.lim,
                      X=X, W=W,W_c = W_c, initU = optimal.u, option = option, control= list('trace'=1))
        #test.brent <- optim(par = init.lambda, fn =  GetVBICOptim, method = "Brent", lower = min(lambdas)*0.1, upper = max(lambdas),
        #              X=X, W=W, initU = optimal.u, option = option,bic.method = 1, ev = 1, control= list('trace'=1))
        v.fits <- FitVs(X,W,W_c, optimal.u,c(test$par),option, weighted = TRUE)
        bic.list.v <- c(test$value)
        lambdas <- c(test$par)
    }else
    {
      v.fits <- FitVs(X,W,W_c, optimal.u,lambdas,option, weighted = TRUE)
      bic.list.v <- unlist(v.fits$BIC[,bic.type])
      #Pick the best choice from here, using threshold.
    }
    optimal.iter.dat <- selectOptimalInstance(v.fits, bic.list.v, lambdas)
    optimal.v <- optimal.iter.dat$m
    rec.dat$lambda.s <- c(rec.dat$lambda.s,optimal.iter.dat$p)

    #PercentVarEx(as.matrix(X)*as.matrix(W), v = optimal.v)
    #message("Updating K, iter ", i)
    rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v))


    #Record new data
    rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
    rec.dat$lambdas <- UpdateAndCheckSparsityParam(rec.dat$lambdas, lambdas, errorReport = TRUE)
    rec.dat$bic.l <- c(rec.dat$bic.l,bic.list.v)
    if(ncol(optimal.u) == ncol(optimal.v))
    {
      rec.dat$X_hat <- c(rec.dat$X_hat, norm(optimal.u %*% t(optimal.v), type = 'F'))
    }else
    {
      rec.dat$X_hat <- c(rec.dat$X_hat, NA)
    }

    #update the parameters for U based on the new V
    u.sparsity <- DefineSparsitySpace(X,W,W_c, optimal.v, "U", option,reg.elements=reg.elements) #Here in case we hit the max.

    if(use.optim)
    {
      alphas <- SelectCoarseSparsityParamsGlobal(u.sparsity, n.points = 15)
    }else
    {
      alphas <- ProposeNewSparsityParams(bic.list.u, alphas, (u.sparsity), i, n.points = 7, no.drop = no.drop.live)
    }

    if(i > min.iter){

      NOT.CONVERGED <- !checkConvergenceBICSearch(i, rec.dat, conv.mode = CONV.MODE) #returns true if convergence is reached

      message("ongoing objective")

    }
    #record the ongoing optimization information.
    jk <- NA
    if(ncol(optimal.u) == ncol(optimal.v))
    {
      jk <- compute_obj(X, W,W_c, optimal.u, optimal.v, option, decomp = TRUE, loglik = TRUE)
    }

    rec.dat$sparsity.obj[[i]] <- jk
    if(ncol(optimal.v) == 1 & option$fixed_ubiq)
    {
      message("Only down to 1 factor, best be stopping now.")
      #align
      aligned <- AlignFactorMatrices(X,W,optimal.u,optimal.v); optimal.u <- aligned$U; optimal.v<- aligned$V
      if(i < min.iter)
      {
        message("Zeroed-out the results very quickly, this suggests some kind of instability...")
      }
      NOT.CONVERGED <- FALSE
    }
    #align
    aligned <- AlignFactorMatrices(X,W,optimal.u,optimal.v); optimal.u <- aligned$U; optimal.v<- aligned$V
    #Check convergence
    #Possible that the tracking of alternative converngcnes is huring our mmeory..
    conv.options[[i]] <- trackAlternativeConvergences(X,W,W_c, optimal.u,optimal.v, i, rec.dat, u.fits$sparsity_space, v.fits$sparsity_space, option)
    #Will this break something?
    i = i+1
  }
  final.index <- i-1
  if(CONV.MODE == "BIC.change")
  {
    #we want to release the previous iteration, not the next

    rec.dat$Vs[[i]] <- optimal.v
    final.index <- i-2
    optimal.v <- rec.dat$Vs[[final.index]]
  }
  #TODO: add in drops for pve
  #This returns all the data from the last iteration
  #TODO: clean this up. This is very confusing.
  return(list("optimal.v" = optimal.v,"resid.var" = u.fits$resid_var,
              "rec.dat" = rec.dat, "lambda"=rec.dat$lambda.s[final.index], "alpha"=rec.dat$alpha.s[final.index], "options" = option,
              "K"= ncol(optimal.v), "alpha.path" = rec.dat$alpha.s, "lambda.path" = rec.dat$lambda.s, "optimal.u" = optimal.u, "convergence.options" = conv.options))
}

#Convergence criteria for the BIC ssearch
#Converges when K is unchanging from one run to the next, and the percentage size change in the alpha/lambda paramters is less than 5%
#Might consider making this more generous- if it stays on the same log scale, then that is probably good enough....
checkConvergenceBICSearch <- function(index, record.data, conv.perc.thresh = 0.05, conv.mode = "SEPARATE")
{
  if(index > 10)
  {
    message("Late stage convergence, terminate soon....")
  }
  conv.K <- record.data$Ks[[index]] == record.data$Ks[[index-1]]
  if(conv.mode == "COMB.SUM")
  {
    message("Need to test and debug still, but okay")
    EPSILON = 1e-4
   sum.prev = record.data$alpha.s[[index-1]] + record.data$lambda.s[[index-1]]
   sum.curr = record.data$alpha.s[[index]] + record.data$lambda.s[[index]]

   return(  conv.K & (sum.prev - sum.curr) < EPSILON)

  }else if(conv.mode == "BIC.change")
    {
    #Goes until the BIC stops dropped, or the change crosses 0.
    #I think maybe you need a better condition.. if you stop too soon you miss a nice minimum somewhere else.
    sum.bic.score <- record.data$bic.a + record.data$bic.l
    delta.bic.score <- sapply(2:length(sum.bic.score), function(i) sum.bic.score[i-1] > sum.bic.score[i])
    li <- length(delta.bic.score)
    #if(any((!delta.bic.score)))
    if(!delta.bic.score[li]) #The last element isn't decreasing
    {
      #CAse 1- we've only been increasing so far, give it some more time
      if(all(!delta.bic.score)[-li]) #all have been increasing
      {
        message("BIC only increasing so far, give it some time...")
        return(FALSE)
      }
      selected.index <- min(which(!delta.bic.score))
      if(selected.index != (index - 1))
      {
        message("behavior not as expected, second to last index isn't the one declining..")
        message("This would happen if it started increasing early on and didn't swa[]")
        print(sum.bic.score)
        print(selected.index)
      }
      return(TRUE)
    }else
    {
      return(FALSE)
    }


  }  else if(conv.mode == "mat.change")
  {
    if(norm(record.data$Vs[[index]], "f") -  norm(record.data$Vs[[index-1]], "f") < 1e-4){return(TRUE)}
    return(FALSE)
  }
  else
  {
    queries <- c(conv.K,
                 abs(record.data$alpha.s[[index]] - record.data$alpha.s[[index-1]])/record.data$alpha.s[[index-1]] < conv.perc.thresh,
                 abs(record.data$lambda.s[[index]] - record.data$lambda.s[[index-1]])/record.data$lambda.s[[index-1]] < conv.perc.thresh)
    return(all(queries))
  }




}

getBICOptimalIndex <- function(record.data)
{
  sum.bic.score <- record.data$bic.a + record.data$bic.l
  delta.bic.score <- sapply(2:length(sum.bic.score), function(i) sum.bic.score[i-1] > sum.bic.score[i])
  return(which.min(sum.bic.score ))
  #TODO: verify the index here is the same as the index we'd get elsewhere.

}
#TODO: try with random, and with non-random.
#gwasML_ALS_Routine(opath, option, X, W, bic.dat$optimal.v)
# gwasML_ALS_Routine(opath, option, X, W, bic.dat$optimal.v)
gwasML_ALS_Routine <- function(option, X, W,W_c, optimal.init, maxK=0, opath = "", no.prune = FALSE,reg.elements=NULL,...)
{

  message(paste0("Starting at k:", option$K))
#if(option$u_init != "")
#{
#  reg.run <- Update_FL(X, W, option, preU = optimal.init)
#}
option$debug <- TRUE
reg.run <- Update_FL(X, W, W_c, option, preV = optimal.init, reg.elements =reg.elements,... )
if(no.prune){return(reg.run)}
if(maxK != 0)
{
  if(ncol(reg.run$V)> maxK)
  {
    message("Parsing down to desired number of factors")
    message("More columns exist than specified....")
  }

  or <- sort(reg.run$PVE, decreasing =TRUE, index.return= TRUE)
  if(ncol(reg.run$V) < maxK)
  {
    message("Resulted in fewer than desired columns. Sorry.")
    maxK <- ncol(reg.run$V)
  }
}else
{
  maxK <- ncol(reg.run$V)
}
  reg.run <- PruneNumberOfFactors(X,W,W_c,reg.run,option$Kmin, maxK, option) %>% OrderEverythingByPVE(.,X,W)


save(reg.run, file = paste0(option$out,opath, "_gwasMF_iter.Rdata" ))
#print(paste0(option$out,opath, "_gwasMF_iter.Rdata" ))

if(option$plots)
{
  o <- data.frame("rownames" = colnames(X), reg.run$V)
  write.table(o, file =  paste0(option$out,opath, ".factors.txt"), quote= FALSE, row.names = FALSE)
  source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/plot_functions.R")
  title <- paste0("A=",reg.run$autofit_alpha[1], "Lambda=",reg.run$autofit_lambda[1])
  plotFactors(apply(reg.run$V, 2, function(x) x/norm(x, "2")), trait_names = o$rownames, title = title)
  ggsave(paste0(option$out,opath, ".factors.png"))
}
#DEBUGGING: objective differs:
return(reg.run)
}

#runFullPipeClean(args$prefix,args, gwasmfiter =args$bic_adj)
runFullPipeClean <- function(args, gwasmfiter =5, save.pipe = FALSE,rep.run = FALSE,...)
{
  opath = "BIC_"
  message("Current setting is to use U based sparsity each time...")
  option <- readInSettings(args)
  option$swap <- FALSE
  option$alpha1 <- 1e-10
  option$lambda1 <- 1e-10
  output <- args$output
  bii = 0
  #Read in the hyperparameters to explore
  hp <- readInParamterSpace(args)
  input.dat <- readInData(args)
  X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names; W_c <- input.dat$W_c
  if(option$K == 0)
  {
    message('Iniitializing X to the max -1')
    option$K <- ncol(X)-1
  }
  #Run the bic thing...
  option$svd_init <- TRUE
  option$fixed_ubiqs <- TRUE
  if(rep.run)
  {
    set.seed(22)
    option$svd_init <- FALSE
    bii = 4
    bic.dat2 <- getBICMatrices(opath,option,X,W,W_c,all_ids, names, burn.in.iter = bii, ...)
    bic.dat3 <- getBICMatrices(opath,option,X,W,W_c,all_ids, names, burn.in.iter = bii, ...)
    bic.dat4 <- getBICMatrices(opath,option,X,W,W_c,all_ids, names, burn.in.iter = bii, ...)
    bic.dat5 <- getBICMatrices(opath,option,X,W,W_c,all_ids, names, burn.in.iter = bii, ...)
      if(save.pipe) {
        save(bic.dat2,bic.dat3,bic.dat4,bic.dat5, file = paste0(option$out,"BIC_iter.MULTIPLES.Rdata" ))
      }
    }


  bic.dat <- getBICMatrices(opath,option,X,W,W_c, all_ids, names, burn.in.iter = bii, ...)

  if(save.pipe) {
  save(bic.dat, file = paste0(option$out,opath, "BIC_iter.Rdata" ))
  }
  option <- bic.dat$options
  option$K <- bic.dat$K
  option$alpha1 <- bic.dat$alpha
  option$lambda1 <- bic.dat$lambda

  ret <- gwasML_ALS_Routine(option, X, W, W_c, bic.dat$optimal.v, maxK=bic.dat$K) #I like this better

  ret[["snp.ids"]] <- all_ids
  ret[["trait.names"]] <- names
  list("als" = ret, "bic.search"=bic.dat)
}




#' Wrapper to run
#'
#' @param args List of all arguments needed to run the pipeline, including path to input data files
#' @param alpha \eqn{\alpha} sparsity parameter for U matrix
#' @param lambda \eqn{\lambda} sparsity parameter for V matrix (numeric)
#' @param opath (optional) prefix for output files (string)
#'
#' @return a list containing all return gwasMF items
#' @export
runStdPipeClean <- function(args,alpha,lambda, opath = "", initV = NULL)
{
  option <- readInSettings(args)
  option$swap <- FALSE
  output <- args$output
  #Read in the hyperparameters to explore
  hp <- readInParamterSpace(args)
  input.dat <- readInData(args)
  X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names; W_c <- input.dat$W_c
  if(option$K == 0)
  {
    message('Iniitializing X to the max')
    option$K <- ncol(X)
  }

  option$alpha1 <- alpha
  option$lambda1 <- lambda
  ret <-gwasML_ALS_Routine(option, X, W, W_c, initV, maxK = option$K)

  ret[["snp.ids"]] <- all_ids
  ret[["trait.names"]] <- names
  ret
}


#' Wrapper to programmatically call the BIC-based factorization
#'
#' @param X data matrix to factorize, SNPs x studies
#' @param W uncertainty weights (1/SE), same dimensions as X
#' @param snp.ids the names of the SNPs in X, in order
#' @param trait.names the names of the studies in X, in order
#' @param K the target number of K to yield.
#' @param gwasmfiter how many burn-in iterations to run, only used if rep.run is TRUE
#' @param rep.run if should proceed with multiple random initializations. Otherwise, initalize with SVD
#' @param bic.var How to calculate the variance on BIC term. Deprecated as of Feb 21, 2023
#' @param use.init.k To initialize to given K or to full
#' @param init.mat Which matrix to initialize with, V or U. Default is V
#' @param is.sim Specify if this run is a simulation. Default is false
#' @param save.path save results along the way as output (optional)
#'
#' @return
#' @export
#' #res <- gleaner(true.betas,1/true.ses, snps, colnames(true.ses), K=5,C=NULL)
gleaner <- function(X,W, snp.ids, trait.names, C = NULL, K=0, gwasmfiter =5, rep.run = FALSE, covar_se=NULL,
                    bic.var= "sklearn", use.init.k = FALSE, init.mat = "V", is.sim = FALSE,
                    save.path = "", scale.mats = FALSE, regression.method = "glmnet", shrinkWL=-1,...)
{
  opath = ""
  if(is.null(C))
  {
    C = diag(ncol(X))
  }


  d <- initializeGwasMF(X,W,C, snp.ids, trait.names, K=ifelse(use.init.k, K, 0),
                        init.mat=init.mat, covar_shrinkage=shrinkWL,covar_se=covar_se,...) #Either use specified, or prune down as we
  option <- d$options; args <- d$args; hp <- d$hp; all_ids <- d$all_ids; names <- d$namesl; W_c <- d$W_c
  print(W_c)
  #pick up here.
  option$block_covar <- 0.2
  #option$WLgamma <- 0.6
  if(is.sim)
  {
    option$Kmin <- K
  }
  option$scale <- scale.mats
  message("scaling set to: ", option$scale)
  option$svd_init <- TRUE; args$svd_init <- TRUE
  K.cap <- K
  option$bic.var <- bic.var
  print(paste("bic var setting:", bic.var))
  #Do the bic learning...
  use.optim <- TRUE
  option$regression_method = regression.method
  if(option$regression_method == "glmnet")
  {
    bic.dat <- getBICMatricesGLMNET(opath,option,X,W,W_c, all_ids, names)
  }else
  {
    bic.dat <- getBICMatrices(opath,option,X,W,W_c, all_ids, names)
  }
  if(is.sim)
  {
    save(bic.dat, file = paste0(save.path, "bic.RData"))
  }

  if(rep.run)
  {
    option$svd_init <- FALSE
    bic.dat2 <- getBICMatrices(opath,option,X,W,W_c, all_ids, names)
    bic.dat3 <- getBICMatrices(opath,option,X,W,W_c,all_ids, names)
    bic.dat4 <- getBICMatrices(opath,option,X,W,W_c,all_ids, names)
    bic.dat5 <- getBICMatrices(opath,option,X,W,W_c,all_ids, names)
  }

  if(is.na(bic.dat$K) | is.na(bic.dat$alpha))
  {
    message("No signal detected under current settings. Program will end")
     return(list("V" = matrix(0, nrow=ncol(X), ncol = 1), "U" = matrix(0, nrow=nrow(X), ncol = 1), "initK" = NA, "K" =0, "obj" = c(NA), "obj_change" = c(),
                     "V_sparsities" = c(), "U_sparsities" = c(), "autofit_lambda" = c(), "autofit_alpha"=c(), "mse" = c(),
                     "V_change" = c(), "U_change" = c(), "decomp_obj" ="", "model.loglik" = c(), "Vs"=list(), "Us"=list(), "pve"=c(0)))
  }
  option <- bic.dat$options
  option$V <- TRUE
  option$alpha1 <- bic.dat$alpha #this is used to regularize U
  option$lambda1 <- bic.dat$lambda #This is used to regularize V
  option$K <- bic.dat$K
  if(rep.run)
  {
    option$K <-  ncol(X)-1
    #message("resetting to full k..")
    alphas.many <- c(bic.dat$alpha, bic.dat2$alpha, bic.dat3$alpha, bic.dat4$alpha, bic.dat5$alpha)
    lambdas.many <- c(bic.dat$lambda, bic.dat2$lambda, bic.dat3$lambda, bic.dat4$lambda, bic.dat5$lambda)
    all.runs <- list()
      for(i in 1:length(alphas.many))
      {
        print(i)
        ret.list <- list()
        for(j in 1:5)
        {
          print(j)
          option <- bic.dat$options
          # option$K <- bic.dat$K
          option$K <-  ncol(bic.dat2$optimal.v)
          #message("resetting to full k..")
          option$V <- TRUE
          option$alpha1 <- alphas.many[i] #Best from the previous. Seems to have converd..
          option$lambda1 <- lambdas.many[i]
          ret.list[[j]] <- gwasML_ALS_Routine(option, X, W,W_c, NULL, maxK=K.cap)
        }
        all.runs[[i]] <- ret.list
      }

    #Select the one that is most stable
    source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/cophenetic_calc.R")
    #move into simple list
    l.sub <- lapply(all.runs, function(ll) lapply(ll, function(x) x$V))
    cophs <- sapply(l.sub, function(x) copheneticAssessment(x))
    print(cophs)
    best <- which.max(cophs)
    #run the max again
    ret <- all.runs[[best]][[1]]
  } else
  {
    #test 1: 6 in, drops 1 CHECK
    #test 2: 2 in, drops none CHECK
    #test3: 8 in, drops 3
    ret <- gwasML_ALS_Routine(option, X, W, W_c, bic.dat$optimal.v, maxK=K.cap)
    if(is.sim)
    {
      save(ret, file = paste0(save.path, "global.fit.als.RData"))
    }
    #Try with 1st params
    #option$alpha1 <- bic.dat$alpha.path[1] #this is used to regularize U
    #option$lambda1 <- bic.dat$lambda.path[1] #This is used to regularize V
  }
   #randomly initialize here...
  #ret.std <- gwasML_ALS_Routine(opath, option, X, W, NULL, maxK=initk)
  #ret2 <- gwasML_ALS_Routine(opath, option, X, W, NULL, maxK=initk) #randomly initialize here...
  #ret3 <- gwasML_ALS_Routine(opath, option, X, W, NULL, maxK=initk) #randomly initialize here...
  ret
}

#####
#Version optimized for glmnet
#bic.dat <- getBICMatricesGLMNET(opath,option,X,W,W_c, all_ids, names)
getBICMatricesGLMNET <- function(opath,option,X,W,W_c, all_ids, names, min.iter = 2, max.iter = 100, burn.in.iter = 0,
                                 init.mat = NULL, rg_ref=NULL,reg.elements=NULL)
{
  W_ld = NULL
  conv.options <- list()
  #Data to track
  rec.dat <- list("alphas"=c(), "lambdas"=c(), "bic.a" = c(), "bic.l"=c(), "obj"=c(),
                  "v.sparsity" = c(), "u.sparsity"=c(), "iter"=c(), "sd.sparsity.u" = c(), "sd.sparsity.v" = c(),
                  "alpha.s" = c(), "lambda.s" = c(), "Ks" = c(), "Vs" = list(), "X_hat" = c(), "sparsity.obj" <- list())
  #If we get columns with NA at this stage, we want to reset and drop those columns at the beginning.
  #Currently, just using NULL for W_ld
  #In the case we don't want any sparsity, just terminate early.
  #print(pryr::mem_used())
  print("Define sparsity space")
  burn.in.sparsity <- DefineSparsitySpaceInit(X, W, W_c, NULL, option, burn.in = burn.in.iter, rg_ref = rg_ref,reg.elements=reg.elements) #If this finds one with NA, cut them out here, and reset K; we want to check pve here too.
  #print(paste0("Current memory usage:", pryr::mem_used()))
  #If we don't want to get the BIC
  if(option$bic.var == "NONE")
    {
      return(list("optimal.v" = burn.in.sparsity$V_burn,"resid.var" = NA,
                  "rec.dat" = rec.dat, "lambda"=1e-20, "alpha"=1e-20, "options" = option,
                  "K"= burn.in.sparsity$new.k, "alpha.path" = NA, "lambda.path" = NA, "optimal.u" = NA, "convergence.options" = NA))
      #We don't need U passed in unless it is request to innitialize U, but not considering for th emoment.
    }


  if(option$u_init != "")
  {
    optimal.u <- burn.in.sparsity$U_burn
    option$K <- burn.in.sparsity$new.k
    message("initializing with U")
  }else if(!is.null(init.mat))
  {
    optimal.v <- init.mat
    option$K <- ncol(optimal.v)
  }else
  {
    optimal.v <- burn.in.sparsity$V_burn
    option$K <- ncol(optimal.v)
  }

  #things to record

  INCR.LIMIT=3
  NOT.CONVERGED <- TRUE; i = 1

  CONV.MODE <- option$param_conv_criteria
  message("Convergence mode: ", CONV.MODE)
  option$regression_method = "glmnet"
  while(i < max.iter & NOT.CONVERGED){
    if(option$u_init != "")
    {
      v.fits <- FitVWrapper(X,W,W_c, optimal.u,option,reg.elements=reg.elements)
      rec.dat$lambda.s <- c(rec.dat$lambda.s,v.fits$lambda.sel)
      optimal.v <- v.fits$V
      optimal.v <- DropEmptyColumns(optimal.v)
      rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v)) #this is based on the previous model?, we haven't pruned out empty columns yet?

      #Record new data
      rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
      rec.dat$bic.l <- c(rec.dat$bic.l,unlist(v.fits$BIC[,bic.type]))
      #update the parameters for U based on the new V
      u.sparsity <- DefineSparsitySpace(X,W,W_c, optimal.v, "U", option,reg.elements=reg.elements) #Here in case we hit the max.
      alphas <- SelectCoarseSparsityParamsGlobal(u.sparsity, n.points = 15)
    }
    #now fit U:
    rec.dat$Vs[[i]] <- optimal.v
    bic.var = option$bic.var
    option$alpha1 <- NA # Important
    option$lambda1 <- NA #important.
    print(paste0("Fitting U, iteration ", i))
    u.fit <- FitUWrapper(X, W, W_c, optimal.v, option,reg.elements=reg.elements)
    #print(pryr::mem_used())

    rec.dat$alphas <- UpdateAndCheckSparsityParam(rec.dat$alphas, u.fit$alpha.sel, errorReport = TRUE)
    rec.dat$bic.a <- c(rec.dat$bic.a,u.fit$bic)
    #Pick the best choice from here, using threshold.
    optimal.u <- u.fit$U
    rec.dat$alpha.s <- c(rec.dat$alpha.s,u.fit$alpha.sel)
    rec.dat$U_sparsities = c(rec.dat$U_sparsities, matrixSparsity(optimal.u, ncol(X)));
    #TODO: drop empty columns here, now.
    optimal.u <- DropEmptyColumns(optimal.u)
    #Do a check- is it empty?
    if(CheckUEmpty(optimal.u)) {
      message("U has no signal; ending model selection phase");
      return(list("optimal.v" = matrix(0, ncol = ncol(optimal.u), nrow = ncol(optimal.v)),"resid.var" = u.fit$resid_var,
                  "rec.dat" = rec.dat, "lambda"=rec.dat$lambda.s[i], "alpha"=rec.dat$alpha.s[i], "options" = option,
                  "K"= NA, "alpha.path" = rec.dat$alpha.s, "lambda.path" = rec.dat$lambda.s, "optimal.u" = optimal.u, "convergence.options" = conv.options))
      }


    #fit V now
    print(paste0("Fitting V, iteration ", i))
    v.fits <- FitVWrapper(X,W,W_c, optimal.u,option,reg.elements=reg.elements)
    #print(pryr::mem_used())
    rec.dat$lambda.s <- c(rec.dat$lambda.s,v.fits$lambda.sel)
    rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v)) #this is based on the previous model?, we haven't pruned out empty columns yet?
    optimal.v <- v.fits$V
    optimal.v <- DropEmptyColumns(optimal.v)
    #Record new data
    rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
    rec.dat$lambdas <- UpdateAndCheckSparsityParam(rec.dat$lambdas, v.fits$lambda.sel, errorReport = TRUE)
    rec.dat$bic.l <- c(rec.dat$bic.l,v.fits$bic)
    if(CheckVEmpty(optimal.v)) {
      message("V has no signal; ending model selection phase");
      return( list("optimal.v" = optimal.v, nrow = ncol(optimal.v),"resid.var" = u.fit$resid_var,
                "rec.dat" = rec.dat, "lambda"=rec.dat$lambda.s[i], "alpha"=rec.dat$alpha.s[i], "options" = option,
                "K"= NA, "alpha.path" = rec.dat$alpha.s, "lambda.path" = rec.dat$lambda.s,
           "optimal.u" =  matrix(0, ncol = ncol(optimal.v), nrow = ncol(optimal.u)), "convergence.options" = conv.options))
    }

    if(ncol(optimal.u) == ncol(optimal.v))
    {
      rec.dat$X_hat <- c(rec.dat$X_hat, norm(optimal.u %*% t(optimal.v), type = 'F'))
    }else
    {
      rec.dat$X_hat <- c(rec.dat$X_hat, NA)
    }

    #update the parameters for U based on the new V
    u.sparsity <- DefineSparsitySpace(X,W,W_c, optimal.v, "U", option,reg.elements=reg.elements) #Here in case we hit the max.
    if(i > min.iter){

      NOT.CONVERGED <- !checkConvergenceBICSearch(i, rec.dat, conv.mode = CONV.MODE,conv.perc.thresh = 0.05) #returns true if convergence is reached

      message("ongoing objective")

    }
    #record the ongoing optimization information.
    #NOTHING BELOW changed, can use same as othr versions of code.
    jk <- NA
    if(ncol(optimal.u) == ncol(optimal.v))
    {
      jk <- compute_obj(X, W,W_c, optimal.u, optimal.v, option, decomp = TRUE, loglik = TRUE) #this is giving issues.
    }

    rec.dat$sparsity.obj[[i]] <- jk

    aligned <- AlignFactorMatrices(X,W,optimal.u,optimal.v); optimal.u <- aligned$U; optimal.v<- aligned$V
    #Check convergence
    conv.options[[i]] <- trackAlternativeConvergences(X,W,W_c, optimal.u,optimal.v, i, rec.dat, u.fit$sparsity_space, v.fits$sparsity_space, option)
    #Will this break something?
    i = i+1
  }
  final.index <- i-1
  if(CONV.MODE == "BIC.change")
  {
    #we want to release the previous iteration, not the next
    rec.dat$Vs[[i]] <- optimal.v
    final.index <- i-2 #Why are we doing this?
    #final.index <- getBICOptimalIndex(rec.dat)
    optimal.v <- rec.dat$Vs[[final.index]]
  }

  if(ncol(optimal.u) != ncol(optimal.v))
  {
    message("Warning: final matrices weren't aligned. Aligning now. Note to developer: check this logical condition.")
    aligned <- AlignFactorMatrices(X,W,optimal.u,optimal.v); optimal.u <- aligned$U; optimal.v<- aligned$V
  }
  #This returns all the data from the last iteration
  #TODO: clean this up. This is very confusing.
  return(list("optimal.v" = optimal.v,"resid.var" = u.fit$resid_var,
              "rec.dat" = rec.dat, "lambda"=rec.dat$lambda.s[final.index], "alpha"=rec.dat$alpha.s[final.index], "options" = option,
              "K"= ncol(optimal.v), "alpha.path" = rec.dat$alpha.s, "lambda.path" = rec.dat$lambda.s, "optimal.u" = optimal.u, "convergence.options" = conv.options))
}

#OLD version. Above is for cleanup.
getBICMatricesGLMNET_DEPRECATED <- function(opath,option,X,W,W_c, all_ids, names, min.iter = 2, max.iter = 100, burn.in.iter = 0,
                                 init.mat = NULL, rg_ref=NULL,reg.elements=NULL)
{
  bic.type = 4
  W_ld = NULL
  conv.options <- list()
  #Data to track
  rec.dat <- list("alphas"=c(), "lambdas"=c(), "bic.a" = c(), "bic.l"=c(), "obj"=c(),
                  "v.sparsity" = c(), "u.sparsity"=c(), "iter"=c(), "sd.sparsity.u" = c(), "sd.sparsity.v" = c(),
                  "alpha.s" = c(), "lambda.s" = c(), "Ks" = c(), "Vs" = list(), "X_hat" = c(), "sparsity.obj" <- list())
  #If we get columns with NA at this stage, we want to reset and drop those columns at the beginning.
  #Currently, just using NULL for W_ld
  #In the case we don't want any sparsity, just terminate early.
  #print(pryr::mem_used())
  print("Define sparsity space")
  burn.in.sparsity <- DefineSparsitySpaceInit(X, W, W_c, NULL, option, burn.in = burn.in.iter, rg_ref = rg_ref,reg.elements=reg.elements) #If this finds one with NA, cut them out here, and reset K; we want to check pve here too.
  #print(paste0("Current memory usage:", pryr::mem_used()))
  if(option$bic.var == "NONE")
  {
    return(list("optimal.v" = burn.in.sparsity$V_burn,"resid.var" = NA,
                "rec.dat" = rec.dat, "lambda"=1e-20, "alpha"=1e-20, "options" = option,
                "K"= burn.in.sparsity$new.k, "alpha.path" = NA, "lambda.path" = NA, "optimal.u" = NA, "convergence.options" = NA))
    #We don't need U passed in unless it is request to innitialize U, but not considering for th emoment.
  }


  if(option$u_init != "")
  {
    optimal.u <- burn.in.sparsity$U_burn
    option$K <- burn.in.sparsity$new.k
    message("initializing with U")
  }else if(!is.null(init.mat))
  {
    optimal.v <- init.mat
    option$K <- ncol(optimal.v)
  }else
  {
    optimal.v <- burn.in.sparsity$V_burn
    #Doesn't appear the order makes a big difference.
    #u.sparsity <- DefineSparsitySpace(X,W,as.matrix(optimal.v),"U", option)
    option$K <- ncol(optimal.v)
  }
  consider.params <- SelectCoarseSparsityParams(burn.in.sparsity, burn.in.iter, n.points = 15)

  #things to record

  INCR.LIMIT=3
  NOT.CONVERGED <- TRUE; i = 1
  #CONV.MODE = "any"
  CONV.MODE <- option$param_conv_criteria
  #CONV.MODE = "BIC.change"
  message("Convergence mode: ", CONV.MODE)
  option$regression_method = "glmnet"
  while(i < max.iter & NOT.CONVERGED){
    #TEST version
    #if(i == 1 & args$WLgamma == "MLE")
    if(FALSE)
    {
      message("Using the iterative update version to estimate gamma. Not recommended.")
      gamma.list <- c(option$WLgamma)
      option$alpha1 <- 1e-20
      option$lambda1 <- 1e-20
      for(j in 1:5)
      {
        u.fit <- FitUWrapper(X, W, W_c, optimal.v, option,reg.elements=reg.elements)
        #Fit V
        v.fits <- FitVWrapper(X,W,W_c, u.fit$U,option,reg.elements=reg.elements)
        #Update Wc
        option$WLgamma <- updateGammaShrinkage(X,W,option$C,u.fit$U,v.fits$V)
        gamma.list <- c(gamma.list, option$WLgamma)
        #apply shrinkage
        new.C <- linearShrinkLWSimple(option$C, option$WLgamma)
        #Update Wc
        W_c <- buildWhiteningMatrix(new.C, dim, blockify = -1)$W_c

        #Repat a few times
      }

    }

    if(option$u_init != "")
    {
      v.fits <- FitVWrapper(X,W,W_c, optimal.u,option,reg.elements=reg.elements)
      rec.dat$lambda.s <- c(rec.dat$lambda.s,v.fits$lambda.sel)
      optimal.v <- v.fits$V
      optimal.v <- DropEmptyColumns(optimal.v)
      rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v)) #this is based on the previous model?, we haven't pruned out empty columns yet?


      #Record new data
      rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
      rec.dat$bic.l <- c(rec.dat$bic.l,unlist(v.fits$BIC[,bic.type]))
      #update the parameters for U based on the new V
      u.sparsity <- DefineSparsitySpace(X,W,W_c, optimal.v, "U", option,reg.elements=reg.elements) #Here in case we hit the max.
      alphas <- SelectCoarseSparsityParamsGlobal(u.sparsity, n.points = 15)
    }
    #now fit U:
    rec.dat$Vs[[i]] <- optimal.v
    bic.var = option$bic.var
    option$alpha1 <- NA # Important
    option$lambda1 <- NA #important.
    print(paste0("Fitting U, iteration ", i))
    u.fit <- FitUWrapper(X, W, W_c, optimal.v, option,reg.elements=reg.elements)
    #print(pryr::mem_used())

    rec.dat$alphas <- UpdateAndCheckSparsityParam(rec.dat$alphas, u.fit$alpha.sel, errorReport = TRUE)
    rec.dat$bic.a <- c(rec.dat$bic.a,u.fit$bic)
    #Pick the best choice from here, using threshold.
    optimal.u <- u.fit$U
    rec.dat$alpha.s <- c(rec.dat$alpha.s,u.fit$alpha.sel)
    rec.dat$U_sparsities = c(rec.dat$U_sparsities, matrixSparsity(optimal.u, ncol(X)));
    #TODO: drop empty columns here, now.
    optimal.u <- DropEmptyColumns(optimal.u)
    #Do a check- is it empty?
    if(CheckUEmpty(optimal.u)) {message("U has no signal; ending model selection phase");
      return(list("optimal.v" = matrix(0, ncol = ncol(optimal.u), nrow = ncol(optimal.v)),"resid.var" = u.fit$resid_var,
                  "rec.dat" = rec.dat, "lambda"=rec.dat$lambda.s[i], "alpha"=rec.dat$alpha.s[i], "options" = option,
                  "K"= NA, "alpha.path" = rec.dat$alpha.s, "lambda.path" = rec.dat$lambda.s, "optimal.u" = optimal.u, "convergence.options" = conv.options))
    }


    #fit V now
    print(paste0("Fitting V, iteration ", i))
    v.fits <- FitVWrapper(X,W,W_c, optimal.u,option,reg.elements=reg.elements)
    #print(pryr::mem_used())
    rec.dat$lambda.s <- c(rec.dat$lambda.s,v.fits$lambda.sel)
    rec.dat$Ks <- c(rec.dat$Ks, ncol(optimal.v)) #this is based on the previous model?, we haven't pruned out empty columns yet?
    optimal.v <- v.fits$V
    optimal.v <- DropEmptyColumns(optimal.v)
    #Record new data
    rec.dat$V_sparsities = c(rec.dat$V_sparsities, matrixSparsity(optimal.v, ncol(X)));
    rec.dat$lambdas <- UpdateAndCheckSparsityParam(rec.dat$lambdas, v.fits$lambda.sel, errorReport = TRUE)
    rec.dat$bic.l <- c(rec.dat$bic.l,v.fits$bic)
    if(CheckVEmpty(optimal.v)) {message("V has no signal; ending model selection phase");
      return( list("optimal.v" = optimal.v, nrow = ncol(optimal.v),"resid.var" = u.fit$resid_var,
                   "rec.dat" = rec.dat, "lambda"=rec.dat$lambda.s[i], "alpha"=rec.dat$alpha.s[i], "options" = option,
                   "K"= NA, "alpha.path" = rec.dat$alpha.s, "lambda.path" = rec.dat$lambda.s,
                   "optimal.u" =  matrix(0, ncol = ncol(optimal.v), nrow = ncol(optimal.u)), "convergence.options" = conv.options))
    }

    if(ncol(optimal.u) == ncol(optimal.v))

    {
      rec.dat$X_hat <- c(rec.dat$X_hat, norm(optimal.u %*% t(optimal.v), type = 'F'))
    }else
    {
      rec.dat$X_hat <- c(rec.dat$X_hat, NA)
    }

    #update the parameters for U based on the new V
    u.sparsity <- DefineSparsitySpace(X,W,W_c, optimal.v, "U", option,reg.elements=reg.elements) #Here in case we hit the max.
    if(i > min.iter){

      NOT.CONVERGED <- !checkConvergenceBICSearch(i, rec.dat, conv.mode = CONV.MODE,conv.perc.thresh = 0.05) #returns true if convergence is reached

      message("ongoing objective")

    }
    #record the ongoing optimization information.
    #NOTHING BELOW changed, can use same as othr versions of code.
    jk <- NA
    if(ncol(optimal.u) == ncol(optimal.v))
    {
      jk <- compute_obj(X, W,W_c, optimal.u, optimal.v, option, decomp = TRUE, loglik = TRUE) #this is giving issues.
    }

    rec.dat$sparsity.obj[[i]] <- jk

    aligned <- AlignFactorMatrices(X,W,optimal.u,optimal.v); optimal.u <- aligned$U; optimal.v<- aligned$V
    #Check convergence
    conv.options[[i]] <- trackAlternativeConvergences(X,W,W_c, optimal.u,optimal.v, i, rec.dat, u.fit$sparsity_space, v.fits$sparsity_space, option)
    #Will this break something?
    i = i+1
  }
  final.index <- i-1
  if(CONV.MODE == "BIC.change")
  {
    #we want to release the previous iteration, not the next
    rec.dat$Vs[[i]] <- optimal.v
    final.index <- i-2

    optimal.v <- rec.dat$Vs[[final.index]]
  }
  #u.norms <- sapply(1:length(conv.options), function(i) conv.options[[i]]$U.norm)
  #v.norms <- sapply(1:length(conv.options), function(i) conv.options[[i]]$V.norm)
  #TODO: add in drops for pve
  #This returns all the data from the last iteration
  #TODO: clean this up. This is very confusing.
  return(list("optimal.v" = optimal.v,"resid.var" = u.fit$resid_var,
              "rec.dat" = rec.dat, "lambda"=rec.dat$lambda.s[final.index], "alpha"=rec.dat$alpha.s[final.index], "options" = option,
              "K"= ncol(optimal.v), "alpha.path" = rec.dat$alpha.s, "lambda.path" = rec.dat$lambda.s, "optimal.u" = optimal.u, "convergence.options" = conv.options))
}

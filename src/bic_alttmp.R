#Deviance_based_BIC
fit
#This may not work either.....
BICglm=function(fit){
  #based on code from https://stackoverflow.com/questions/40920051/r-getting-aic-bic-likelihood-from-glmnet
  tLL <- -deviance(fit) # 2*log-likelihood
  k <- fit$df
  n <- nobs(fit)
  AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
  AIC_ <- -tLL+2*k
  
  log(n)*k - tLL
}
 #looks promising
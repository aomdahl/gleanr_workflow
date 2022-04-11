############################################################################################
# Script containing functions to estimate percent variance explained of matrices
############################################################################################
# Ashton Omdahl, Sept 2021
############################################################################################


pacman::p_load(data.table, tidyr, dplyr, reticulate)

#default setting is to always return an estimate of "additive" PVE, not cumulative.
#As in, what is the gain in explained variance by adding in one more factor?

#Wrapper to get all of them at once in a data frame.
calcPVEAll <- function(V,U,X,K,D = NULL)
{
  if(is.null(D))
  {
    D = rep(1, K)
  }
  var_tracker <- data.frame("Factor" = 1:K)
  var_tracker$r2 <- r2Score(V,U,X,K,D)
  var_tracker$var_exp <- expVarScore(V,U,X,K,D)
  var_tracker$pveBySVD <- pveBySVD(V,U,X,K,D)
  #var_tracker$flashR <- approxPVE(V,U,X,K,D)
  var_tracker$flashRWorks <-  approxPVE_init(X, V, U %*% diag(D))
  return(var_tracker)
}

#Calculate a simple r2, using sklearn's "r2_score" function
r2Score <- function(V,U,X,K,D = NULL)
{
  sklran <- import("sklearn.metrics")
  x_py= r_to_py(X)
  xhat_py = r_to_py((matrix(U[,1]) * D[1]) %*% t(matrix(V[,1])))
  r2 <- c(sklran$r2_score(x_py, xhat_py), sapply(2:K, function(i) sklran$r2_score(x_py, r_to_py(U[,1:i] %*% diag(D[1:i]) %*% t((V[,1:i]))))))
   return(c(r2[1], sapply(2:K, function(i) r2[i] - r2[i-1])))
}



## Method 2: Explained Variance Score
#Basically, this drops the assumption of the data being 0-centered (? still not clear on this point). This would look like $1-\frac{Var[\hat{y}-y]}{Var[y]}$
#for reference, see:
#`https://stats.stackexchange.com/questions/210168/what-is-the-difference-between-r2-and-variance-score-in-scikit-learn`
#and `https://scikit-learn.org/stable/modules/generated/sklearn.metrics.explained_variance_score.html`.
expVarScore <- function(V,U,X,K,D = NULL)
{
  sklran <- import("sklearn.metrics")
     x_py= r_to_py(X)
   xhat_py <- r_to_py(matrix(U[,1] * D[1]) %*% t(matrix(V[,1])))
   var_exp <- c(sklran$explained_variance_score(x_py, xhat_py), sapply(2:K, function(i) sklran$explained_variance_score(x_py, r_to_py(U[,1:i] %*% diag(D[1:i]) %*% t((V[,1:i]))))))
   return( c(var_exp[1], sapply(2:K, function(i) var_exp[i] - var_exp[i-1])))
 }




## Method 3: Shen and Huan 2008 method
#This is implemented in PMA for sparse PCA, and seems to assume that our loading matrix is orthornomal. Not sure if this will work
#`https://www-sciencedirect-com.proxy1.library.jhu.edu/science/article/pii/S0047259X07000887?via%3Dihub`
#`https://rdrr.io/cran/scorer/src/R/regression_metrics.R#sym-explained_variance_score`
pveBySVD <-function(V,U,X,K,D = NULL)
{
  ve <- c()
  full.svd <- svd(X)
  for(k in 1:K)
  {
    k.now <- matrix(V[,1:k], ncol = k)
    projection <- X %*% k.now %*% solve(t(k.now) %*% k.now) %*% t(k.now)
    p.svd <- svd(projection)
    ve <- c(ve, sum(p.svd$d^2))
  }
  r <- ve / sum((full.svd$d)^2) 
  return(c(r[1], sapply(2:K, function(i) r[i] - r[i-1])))
}


#Method 4: like the one done by FlashR
#Accounts for the residual error and the amount of scaling needed on the factor.
scaledMats <- function(f,l)
{
  d = sqrt(colSums(f^2) * colSums(l^2))
  scaled.l <- scale(l, scale = sqrt(colSums(l^2)), center = FALSE)
  scaled.f <- scale(f, scale = sqrt(colSums(f^2)), center = FALSE)
  return(list(d = d, l = scaled.l, f = scaled.f))
}

approxPVE <- function(V,U,X,K,D = NULL)
{
  if(is.null(D))
  {
    D = rep(1,K)
  }
  d = scaledMats(V,U)$d
  print(D)
  L <- U %*% diag(D)
  r2 <- (X - L %*% t(V))^2 #this could be better done if we had 2nd moments.
  N = nrow(X)
  tau_hat <- 1/r2
  d^2/sum(d^2 + sum(1/tau_hat))
}


#FlashR method 
#why the 2 different ones? ppooh
approxPVE_init <- function(X, F, L)
{
  d = scaledMats(F,L)$d
  r2 <- (X - L %*% t(F))^2 #this could be better done if we had 2nd moments.
  N = nrow(X)
  #tau_hat <- N/colSums(r2) #assumes column-specific variance, 
  #tau_hatm <- t(replicate(nrow(X), tau_hat))
  tau_hat <- 1/r2
  d^2/sum(d^2 + sum(1/tau_hat))
}


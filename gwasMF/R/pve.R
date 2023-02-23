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
#reading further on this later on, I think it works because we aren' requiring this? should dig into it more...
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

mean.na <- function(vec){
  return(mean(vec[!is.na(vec)]))
}

PercentVarEx <- function(x,v, K=NULL, center= FALSE)
{
  if(all(v == 0) | nrow(v) <= 1)
  {
	message("empty data passed in, crisis averted")
	  return(c(NA))
  }
  PMA_PVE(x,v, K=NULL, center= FALSE)
}

getVE <- function(xfill, vk)
{
  new.ve <- tryCatch({
    xk <- xfill%*%vk%*%solve(t(vk)%*%vk, tol = 1e-20)%*%t(vk)
    svdxk <- svd(xk)
    sum(svdxk$d^2)
  },  error=function(cond) {
    message("caught error calculating PVE")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  }
  )
  new.ve
}

PMA_PVE <- function(x,v, K=NULL, center= FALSE)
{
  #Code taken directly from https://rdrr.io/cran/PMA/src/R/PMD.R
  # Calculate percent variance explained, using definition on page 7 of Shen and Huang (2008) Journal of Multivariate Analysis vol 99: 1015-1034
  #ADdition on 1/10: Account for cases with columns of all 0.
  #should report a PVE of 0
  k.init <- ncol(v)
  drops <- c()
  if(any(colSums(v==0) == nrow(v)))
  {
    message("Warning: rows in current v contain 0 entries. Setting PVE to 0.")
    drops <- which(colSums(v==0) == nrow(v))
    v <- v[,-drops]
  }
  if(is.null(K))
  {
    K=ncol(v)
  }
  if(!is.numeric(K))
  {
	  print("K is not numeric....")
	  return(c(NA))
  }
  v <- matrix(v, ncol=K)
  if(ncol(v) == 1 & nrow(v) == 1)
  {
	  return(c(NA))
  }
  ve <- NULL # variance explained
  xfill <- x
  if(center) xfill <- x-mean.na(x)
  xfill[is.na(x)] <- mean.na(xfill)
  for(k in 1:K){
    vk <- matrix(v[,1:k], ncol=k)
    

    ve <- c(ve, getVE(xfill,vk))

    
  }
  pve <- ve/sum(svd(xfill)$d^2) # proportion of variance explained
  if(K > 1)
  {
    pve <- c(pve[1], sapply(2:K, function(i) pve[i] - pve[i-1]))
  }else
  {
    pve <- c(pve[1])
  }
  #If there were drops, add them back in
  if(length(drops) > 0)
  {
    ret <- rep(0, k.init)
    true.val <- 1
    for(i in 1:k.init)
    {
      if(!(i %in% drops))
      {
       ret[i] <- pve[true.val]
       true.val <- true.val + 1
      }else
      {
        message("Factor ", i, "is dropping..")
      }
    }
    pve <- ret
  }
  #pve
pve
  
}




#debug test
debugTest <- function()
{
  message("A test to see if my function is the same as the one in PMD")
  pca <- svd(scale(X))
  pve <- (pca$d^2/sum(pca$d^2))
  pve.pma <- pveBySVD(pca$v, pca$u, scale(X), K = 55)
  pve.pma.direct <- pveBySVDPMD(scale(X), pca, 55)
  pve.pma.direct.true <- c(pve.pma.direct[1], sapply(2:55, function(i) pve.pma.direct[i] - pve.pma.direct[i-1]))
  plot(pve, pve.pma)
  plot(pve-pve.pma)
  plot(pve.pma,pve.pma.direct.true )
  plot(pve.pma.direct.true - pve.pma)
  message("And so it appears to be!")
}

#
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


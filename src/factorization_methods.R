pacman::p_load(data.table, tidyr, dplyr, magrittr, readr, ggplot2, stringr, penalized, cowplot, parallel, flashr)
#This script has function calls to run the different factorization algorithms.
#It includes:
# * sparsePCA
# * stdPCA
# * PMA, cv on one paramter
# * PMA, cv on 2 parameters
# * flashR (with average factor, point normal, adaptive shrinkage, backfitting)
#Note that for the sake of these simulations, W is believed to be 1 across all.
#This might be something worth considering --> we could make a W based on an assigned sample size for each N and an assigned MAf for each SNP


#F, L and X
#x refers to the estimated reconstruction of X. It varies from method to method, depending on if it pulls out a d, etc.
retDat <- function(F, L, X)
{
  ret_list <- list()
  ret_list$F <- F
  ret_list$L <- L
  ret_list$X_hat <- X
  return(ret_list)
}

#Perform Sparse PCA.
#TODO: figure out how to select the optimal alpha/beta parameters here. Some kind of cross-validation
sparsePCA <- function(X, K)
{
  results <- spca(X, k=K, alpha=0.001, beta=1e-3, center=TRUE, scale=TRUE)
  return(retDat(results$loadings, results$scores))
}

#Code 
stdPCA <- function(X, K)
{
  #results <- prcomp(X, rank = K, scale. = TRUE, center =  TRUE)
  scaled <- scale(X)
  results <- svd(X,nu = K, nv = K)
  return(retDat(results$v, results$v, results$u %*% diag(results$d[1:K]) %*% t(results$v)))
}

#Single version, with one set of cross validation
singlePMA <- function(X, K)
{
  suppressWarnings(library(PMA))
  v.out <- PMD.cv(X, type="standard", sumabss=seq(0.2,1, len=20)) #Unclear to me how to set these numbers. Just trying the largest range that works...
  pmas = PMD(X, sumabs = cv.out$bestsumabs, K = K) ## sumabs: sparsity of v, sumabsu: sparsity of u
  return(retDat(pmas$v, pmas$u, pmas$u %*% diag(pmas$d) %*% t(pmas$v)))
}

#Taken directly from yuan
doublePMA <- function(X, K)
{
  suppressWarnings(library(PMA))
  for(sv in seq(1,sqrt(K))){
    for(su in seq(1,10)){
      cv = tune_cor_PMA(as.matrix(X), rankK, su, sv)
      result = rbind(result, c(sv, su, cv))
    }
  }
  ops = result[which.min(result[,3]),]
  sv = ops[1]
  su = ops[2]
  pmas = PMD(as.matrix(X), K=rankK, sumabs = NULL, sumabsu = su, sumabsv = sv)
  return(retDat(pmas$v, pmas$u, pmas$u %*% diag(pmas$d) %*% t(pmas$v)))
}

#just a standard flashR with backfitting
flashRBackfit <- function(X, K)
{
  suppressWarnings(library(flashr))
  suppressWarnings(library(ashr))
  f_ashr = flash(X, backfit = T, greedy = T, Kmax = K)
  l_hat %*% diag(out$d) %*% t(f_hat)
  return(retDat(f_ashr$ldf$f, f_ashr$ldf$l, f_ashr$ldf$l %*% diag(f_ashr$ldf$d) %*% t(f_ashr$ldf$f) ))
}

#FlashR with adaptive shrinkage and backfitting
flashRashR <- function(X, K)
{
  suppressWarnings(library(flashr))
  suppressWarnings(library(ashr))
  flash_object <- flash_set_data(Y = X)
  f_ashr <- flashr:::flash_greedy_workhorse(flash_object,
                                            ebnm_fn = "ebnm_ash",
                                            verbose_output = "odF", Kmax = K)
  
  f_ashr_backfit = flash_backfit(flash_object,f_ashr, ebnm_fn = "ebnm_ash")
  plotFactors(f_ashr_backfit$ldf$f,trait_names =n, title = "General")
  return(retDat(f_ashr_backfit$ldf$f, f_ashr_backfit$ldf$l))

}

#This one should have a fixed number of factors  at K for sure, uses a PMA based initialization
flashRFixK <- function(X, K)
{
  suppressWarnings(library(flashr))
  suppressWarnings(library(ashr))
  flash_object <- flash_set_data(Y = X)
  f = flash_add_factors_from_data(flash_object,K=K)
  f_bf = flash_backfit(flash_object,f)
  return(retDat(f_bf$ldf$f, f_bf$ldf$l))

}

#add an average factor on the front end along with backfitting.
flashierUbiq <- function(X, K)
{
  suppressWarnings(library(flashier))
  suppressWarnings(library(ashr))
  flashier_mix <- flashier::flash(data = X, greedy.Kmax = K, prior.family = list(prior.nonzero.mode(), prior.normal.scale.mix()), 
                                  verbose.lvl = 1, backfit = TRUE)
  factors <- flashier_mix$loadings.pm[[2]]
  loadings <- flashier_mix$loadings.pm[[1]]
  return(retDat(factors, loadings))

}


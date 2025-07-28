pacman::p_load(data.table, tidyr, dplyr, magrittr, readr, ggplot2, stringr, penalized, cowplot, parallel, flashr, PMA)
if("flashier" %in% rownames(installed.packages()) == FALSE) {devtools::install_github("willwerscheid/flashier", build_vignettes = FALSE); library(flashier)}
#This script has function calls to run the different factorization algorithms.


#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/R/data_decorrelation.R")
library(gleanr)
#It includes:
# * sparsePCA
# * stdPCA
# * PMA, cv on one paramter
# * PMA, cv on 2 parameters
# * flashR (with average factor, point normal, adaptive shrinkage, backfitting)
#Note that for the sake of these simulations, W is believed to be 1 across all.

runSingle <- function(method, X, K, z_path, n_path,se_m = 1, covar = NULL,savepath="./",...){
  Xin = X
  if(is.null(se_m))
  {
    W.mat <- NULL
  }else{
    #print(se_m)
    W.mat <- 1/se_m
  }

  message("now starting: ", method)
  switch(method,
         "PMA2" = doublePMA(Xin, K),

         "flashUbiqFilter" = flashierUbiq(Xin, K, W = W.mat, filter = TRUE),
         "flashUbiq_noscale" = flashierUbiq(Xin, K=K),
         "ashr" = flashRashR(Xin, K, W = W.mat),
         "backfit_noscale" = flashRBackfit(Xin, K, W = NULL),
         "backfit" = flashRBackfit(Xin,K, W = W.mat),
         "FLASH" = flashRBackfit(Xin,K, W = W.mat),
         "backfit_SE" = flashRBackfitScaled(Xin, K, W=W.mat),
         "FLASH_SE" = flashRBackfitScaled(Xin, K, W=W.mat),
         "backfit_SE_kroneker" = flashRBackfitScaled(Xin, K, W=W.mat, var.type.param = c(1,2)),
         "FLASH_SE_kroneker" = flashRBackfitScaled(Xin, K, W=W.mat, var.type.param = c(1,2)),
         "flashColumnVar" = flashRTraitVar(Xin, K),
         "FactorGo" = runFactorGo(Xin, K, z_path, n_path, savepath), #llast 2 args for factorgo only)
         "sSVD" = sparseSVD(Xin, K),
         "SVD_beta"= vanillaPCA(Xin,K), #this case its not scaled beforehand, based on upper script.
         "ssvd" = sparseSVD(Xin, K),
         "sPCA" = sparsePCA(Xin, K),
         "PCA" = vanillaPCA(Xin,K), #Xin is scaled upstream, don't worry about this! vERified
         "SVD" = vanillaPCA(Xin,K), #Xin is scaled upstream, don't worry about this!,
         "PCA_chooseK" = vanillaPCA(Xin,NA),
         "SVD_whiten" = whitenedPCA(Xin, K, C=covar),
         "gwasMF_BIC" = runGWASMFBeta(Xin, W=W.mat, K, C=covar,...),
         "GLEANER" = runGWASMFBeta(Xin, W=W.mat, K, C=covar,regression.method = "glmnet",...),
         "GLEANER_glmnet" = runGWASMFBeta(Xin, W=W.mat, K, C=covar,regression.method = "glmnet",...),
         "GLEANER_glmnet_noCovar" = runGWASMFBeta(Xin, W=W.mat, K, C=NULL,regression.method = "glmnet",...),
         "gwasMF_grid" = runGWASMFGrid(Xin, W.mat, K, C=covar,...),
         "gwasMF_BIC_noCovar" = runGWASMFBeta(Xin, W=W.mat, K, C=NULL,...),
         "GLEANER_noCovar" = runGWASMFBeta(Xin, W=W.mat, K, C=NULL,...)

  )

}


#This might be something worth considering --> we could make a W based on an assigned sample size for each N and an assigned MAf for each SNP


#F, L and X
#x refers to the estimated reconstruction of X. It varies from method to method, depending on if it pulls out a d, etc.
#@param pve: percent variance explained per component
#@param XRef: the actual reconstruction of X used, since it varies
#@param X: the predicted X from the output (Xhat)
#@param U: The factor matrix
#@param V: The loading matrix
retDat <- function(V, U, X, Xref, pve = NULL)
{
  ret_list <- list()
  ret_list$V <- V
  ret_list$U <- U
  ret_list$X_hat <- X
  ret_list$X <- Xref
  if(length(pve) == 0)
  {
    pve <- NULL
  } else if(length(pve) > 1)
  {
    pve <- pve
  }
  else if(length(pve) == 1 & !is.null(pve) & is.na(pve))
  {
    pve <- NULL
  }
  else
  {
    message("weird pve...")
    print(pve)
  }
  ret_list$pve <- pve
  return(ret_list)
}

#Perform standard PCA
vanillaPCA <- function(X,K, covar=NULL, whiten=FALSE)
{
  #Now run it on PCA
  choose.k = FALSE
  if(is.na(K) | K == 0)
  {
    K = ncol(X) - 1
    choose.k = TRUE
  }
  if(!is.null(covar) & whiten)
  {
    X <- adjustMatrixAsNeeded(X, covar)
  }

  pca = RSpectra::svds(as.matrix(X), k = K)
  if(choose.k)
  {
    keep.vals <- which(pca$d > mean(pca$d))
  }else
  {
    keep.vals <- 1:K
  }
  v <- as.matrix(pca$v[,keep.vals])
  u <- as.matrix(pca$u[,keep.vals])
  if(length(keep.vals) == 1)
  {
    d <- pca$d[keep.vals]
  }else
  {
    d <- as.matrix(diag(pca$d[keep.vals]))
  }

  return(retDat(v, u, u %*% d %*% t(v), X, pca$d[keep.vals]^2/sum(pca$d^2)))
}

whitenedPCA <- function(X, K, C)
{
  vanillaPCA(X,K, covar=C, whiten=TRUE)
}




#Perform Sparse PCA.
#TODO: figure out how to select the optimal alpha/beta parameters here. Some kind of cross-validation
sparsePCA <- function(X, K)
{
  suppressWarnings(library(sparsepca))
  results <- spca(X, k=K, alpha=0.001, beta=1e-3, center=FALSE, scale=FALSE) #why don't we want to scale again? was it becasue of the comparison?
  return(retDat(results$loadings, results$scores,results$scores %*% t(results$loadings), X, results$eigenvalues^2/sum(results$eigenvalues^2)))
}
#rapid redo- from yuan
sparseSVD <- function(X, K)
{
    message("ssvd is not compatible with R 4- going to use irlba ssvd")
    suppressWarnings(library(irlba))
	ssvds = suppressMessages(ssvd(as.matrix(X), k=K))
	#suppressWarnings(library(ssvd))
	#    ssvds = ssvd(as.matrix(X), method = 'method', r = K)
    return(retDat(ssvds$v, ssvds$u, ssvds$u %*% diag(diag(ssvds$d)) %*% t(ssvds$v), X, diag(ssvds$d^2/sum(ssvds$d^2))))
}
#Code
stdPCA <- function(X, K)
{
  #results <- prcomp(X, rank = K, scale. = TRUE, center =  TRUE)
  #scaled <- scale(X)
  results <- svd(X,nu = K, nv = K)
  return(retDat(results$v, results$u, results$u %*% diag(results$d[1:K]) %*% t(results$v), X, results$d^2/sum(results$d^2)))
}

#Code
sparcePCA_PMA <- function(X, K)
{
  suppressWarnings(library(PMA))
  cv.out <- SPC.cv(X)

  results <- SPC(X, sumabsv=cv.out$bestsumabsv, K = K)
  #return(retDat(results$v, results$u, results$u %*% diag(results$d) %*% t(results$v), X, results$d^2/sum(results$d^2)))
  return(retDat(results$v, results$u, results$u %*% diag(results$d) %*% t(results$v), X,   results$prop.var.explained))
}


#Single version, with one set of cross validation
singlePMA <- function(X, K)
{
  suppressWarnings(library(PMA))
  cv.out <- PMD.cv(X, type="standard", sumabss=seq(0.2,1, len=20)) #Unclear to me how to set these numbers. Just trying the largest range that works...
  pmas = PMD(X, sumabs = cv.out$bestsumabs, K = K) ## sumabs: sparsity of v, sumabsu: sparsity of u
  #return(retDat(pmas$v, pmas$u, pmas$u %*% diag(pmas$d) %*% t(pmas$v), X, pmas$d^2/sum(pmas$d^2)))
  return(retDat(pmas$v, pmas$u, pmas$u %*% diag(pmas$d) %*% t(pmas$v), X, NA))
}

#Taken directly from yuan
doublePMA <- function(X, K)
{
  suppressWarnings(library(PMA))
  result = NULL
  for(sv in seq(1,sqrt(K))){
    for(su in seq(1,10)){
      cv = suppressMessages(tune_cor_PMA(as.matrix(X), K, su, sv))
      result = rbind(result, c(sv, su, cv))
    }
  }
  ops = result[which.min(result[,3]),]
  sv = ops[1]
  su = ops[2]
  pmas = suppressMessages(PMD(as.matrix(X), K=K, sumabs = NULL, sumabsu = su, sumabsv = sv))
  #return(retDat(pmas$v, pmas$u, pmas$u %*% diag(pmas$d) %*% t(pmas$v), X, pmas$d^2/sum(pmas$d^2)))
  return(retDat(pmas$v, pmas$u, pmas$u %*% diag(pmas$d) %*% t(pmas$v), X, NA))
}

#just a standard flashR with backfitting
flashRBackfit <- function(X, K, W = NULL, filter = FALSE)
{
  #suppressWarnings(library(flashr))
  suppressWarnings(library(ashr))
  #f_ashr = suppressMessages(flashr::flash(X, backfit = TRUE, greedy = TRUE, Kmax = K))
  #return(retDat(f_ashr$ldf$f, f_ashr$ldf$l, f_ashr$fitted_values,X, f_ashr$pve))
  if(!is.null(W))
  {
    X <- X*W
  }else{
    print("Not scaling by W, is this right?")
  }
  f_ashr = flashier::flash(X,  greedy_Kmax = K, backfit = TRUE)
  if(f_ashr$n_factors == 0)
  {
    n = nrow(X)
    m = ncol(X)
    return(retDat(matrix(0,m,K), matrix(0,n,K), matrix(0,n,m),X, 0))
  }
  #return(retDat(f_ashr$F_pm, f_ashr$L_pm, f_ashr$L_pm %*% t(f_ashr$F_pm),X, f_ashr$pve))
  return(retDat(f_ashr$F_pm, f_ashr$L_pm, fitted(f_ashr),X, f_ashr$pve))

}

#The variant with Scaling the right way....
#Updated var.tpye.param to be column-wise on sept 25
flashRBackfitScaled <- function(X, K, W = NULL, filter = FALSE, var.type.param = 0)
{
  fscaled = flashier::flash(X, S=1/W, greedy_Kmax = K, backfit = TRUE,var_type = var.type.param)
  if(fscaled$n_factors == 0)
  {
    n = nrow(X)
    m = ncol(X)
    return(retDat(matrix(0,m,K), matrix(0,n,K), matrix(0,n,m),X, 0))
  }
  return(retDat(fscaled$F_pm, fscaled$L_pm, fitted(fscaled),X, fscaled$pve))
}

##To run factor go
runFactorGo <- function(X, K, zpath, npath, savepath)
{
  message("Running factor go via reticulate")
  reticulate::use_condaenv("/home/aomdahl1/.conda/envs/std")
  #import("factorgo")
  reticulate::source_python("/scratch16/abattle4/ashton/snp_networks/factorGo/FactorGo/src/factorgo/cli_mod.py")
  #In simple empirical test, unscaled is better
  command=paste0(zpath, " ", npath, " -k", K, " -o ", savepath)
  run_main(command)
  #run_main("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs//V6_U6_mafmixed_n100000.high_covar_1block_cont_scaling/fake.z.txt /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs//V6_U6_mafmixed_n100000.high_covar_1block_cont_scaling/fake.N.txt -k 5 -o /scratch16/abattle4/ashton/snp_networks/factorGo/sandbox/test")
  #savepath="/scratch16/abattle4/ashton/snp_networks/factorGo/sandbox/test"
  v <- as.matrix(data.table::fread(paste0(savepath, ".Zm.tsv.gz")))
  u <- as.matrix(data.table::fread(paste0(savepath, ".Wm.tsv.gz")))
  pve <- unlist(data.table::fread(paste0(savepath, ".factor.tsv.gz"))[,3])
  #retDat <- function(V, U, X, Xref, pve = NULL)
  #py_run_file("/scratch16/abattle4/ashton/snp_networks/factorGo/FactorGo/src/factorgo/cli.py -h")

  return(retDat(as.matrix(v), as.matrix(u), u %*% t(v), X,pve=pve))
}

runGWASMFBeta <- function(X, W, K,C,...)
{

  #run.data <- runSingle(meth, effect.matrix, K, se_m=se, covar = c.mat, bic.var = args$bic_var,
  #init.mat = args$init_mat, is.sim = TRUE, save.path = paste0(args$outdir, meth))
  if(is.null(C))
  {
    message("Not using C")
  }
 #source("/scratch16/abattle4/ashton/snp_networks/scratch/testing_gwasMF_code/model_selection/bic_autofit_functions.R")
  #library(gleanr)
  #1/16 dropped these lines, since package working
  #files.source = list.files("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/R/")
  #sapply(paste0("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/R/", files.source), source)

  #Had to update to reflect simulations
 res <- gleanr(X,W, paste0("rs", 1:nrow(X)), paste0("T",1:ncol(X)), K=K,C=C,is.sim=TRUE,...) #recall- if was sim and sklearn, forced it to be GRID
  option <- list();   option$fixed_ubiq <- TRUE; option$alpha1 <- res$autofit_alpha[1]; option$lambda1 <- res$autofit_lambda[1]
  if(nrow(res$V) == 0)
  {
    res$V <- matrix(0, nrow = ncol(X), ncol = 1)
  }
  if(nrow(res$U) == 0)
  {
    res$U <- matrix(0, nrow = nrow(X), ncol = 1)
  }
  #print(compute_obj(X, W, res$U, res$V, option, decomp = TRUE, loglik = 900))
  return(retDat(res$V, res$U, X = res$U %*% t(res$V),Xref = X, res$pve))
} 


runGWASMFGrid <- function(X, W, K,C,...)
  {
  message("Checkpoint2")

    #source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/gwasmf_grid_search.R")
    files.source = list.files("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/R/")
  sapply(paste0("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/R/", files.source), source)
  res <- gwasMFGrid(X,W, paste0("rs", 1:nrow(X)), paste0("T",1:ncol(X)), K=K, Kmax = K, C=C)
    return(retDat(res[[1]]$V, res[[1]]$U, X = res[[1]]$U %*% t(res[[1]]$V),Xref = X, res[[2]]$PVE))
  }

#backfitting, column-based variance
flashRTraitVar <- function(X, K, W = NULL, filter = FALSE)
{
  suppressWarnings(library(flashr))
  suppressWarnings(library(ashr))
  f_ashr = suppressMessages(flashr::flash(X, backfit = TRUE, greedy = TRUE, Kmax = K, var_type = "by_column"))

  return(retDat(f_ashr$ldf$f, f_ashr$ldf$l, f_ashr$fitted_values,X, f_ashr$pve))
}

#FlashR with adaptive shrinkage and backfitting
flashRashR <- function(X, K, W = NULL)
{
  flash_object <- flash_set_data(Y = X, S = W)
  vt = "zero"
  if(is.null(W)){vt = "constant"}
  f_ashr <- suppressMessages(flashr:::flash_greedy_workhorse(flash_object,
                                            ebnm_fn = "ebnm_ash",
                                            verbose_output = "odF", Kmax = K, var_type = vt  ))
  #TODO: allow for column specific variances
  f_ashr_backfit = suppressMessages(flash_backfit(flash_object,f_ashr, ebnm_fn = "ebnm_ash", var_type = vt))
  #plotFactors(f_ashr_backfit$ldf$f,trait_names =n, title = "General")
  return(retDat(f_ashr_backfit$ldf$f, f_ashr_backfit$ldf$l, f_ashr_backfit$fitted_values, X, f_ashr_backfit$pve))

}

#This one should have a fixed number of factors  at K for sure, uses a PMA based initialization
flashRFixK <- function(X, K, W= NULL)
{
  suppressWarnings(library(flashr))
  suppressWarnings(library(ashr))
  flash_object <- flashr::flash_set_data(Y = X, S = W)
  f = flash_add_factors_from_data(flash_object,K=K)
  f_bf = flashr::flash_backfit(flash_object,f)
  return(retDat(f_bf$ldf$f, f_bf$ldf$l, f_bf$fitted_values), X, f_bf$pve)

}

flashCompleteVar <- function(X,K,W)
{
  S = 1/W #standard errors
  X <- effect.matrix
  #how about the urls to the examples here you ninny?

  gtex.bf <- flash(
    X,
    greedy_Kmax = 5,
    backfit = TRUE,
    verbose = 0, var_type = c(1,2)
  )
  #The kroneker variance really mucks things. up...
  #I think these special ones just do worse than vaniall flash honestly.
  S.dim <- ? #specify this
  ebnm.prior <- ebnm::ebnm_point_normal
  var_type = c(1,2)
}

flashCompleteNoVar <- function(X,K,W)
{
  Z = X*W #standard errors
  S.dim <- ? #specify this
    ebnm.prior <- ebnm::ebnm_point_normal
  var_type = c(1,2)
#kroneker variance hasn't been implemented yet


  flash.fix.factors(kset = 1, mode = 2)


}

#add an average factor on the front end along with backfitting.
#If reviewers want it, we can do it.
flashierUbiq <- function(X, K, W = NULL, filter = FALSE,var.type.param=0)
{
  #suppressWarnings(library(flashier))
  #suppressWarnings(library(flashr))
	suppressWarnings(library(ashr))
  suppressWarnings(library(ebnm))
 #was in flashier, dropped to flash
  #flashier_mix <- flashier::flash(data = X, S = 1/W, greedy_Kmax = K, prior.family = list(prior.nonzero.mode(), prior.normal.scale.mix()),
  #                                verbose.lvl = 1, backfit = TRUE,var_type = var.type.param)
  #flashier::flash(data = X, S = 1/W, greedy_Kmax = K, ebnm_fn = list(ebnm_normal, ebnm_point_normal)
  # fscaled = flashier::flash(X, S=1/W, greedy_Kmax = K, backfit = TRUE,var_type = var.type.param)

  # Fit flash model (this is an example)'#For teh list of options, see #ebnm
  library(ebnm)
  library(flashier)

#I think this would do it... first factor on F is dense.
#THIS WORKS- can use this if needed

#Note- it seems like the performance is better if I don't account for S as they specify, but whatever
  #Arguments in favor of doing it- might perform better on simulations, so doesn't look as bad
  #Argument against- the first factor is typically the highest variance, which isn't guaranteed to be the case in my simulations. So it may not help.
fl <- flashier::flash_init(data=X,S=1/W) |>
  flash_greedy(Kmax = 1, ebnm_fn = c(ebnm_point_normal,
                                     ebnm_normal)) |>
  flash_greedy(Kmax = 4, ebnm_fn = c(ebnm_point_normal,
                                     ebnm_point_normal))
  flash_backfit() |>
  flash_nullcheck()
#helpful resources: https://github.com/heyuan7676/ts_eQTLs/blob/master/Extended_Methods/Other_MF_methods/flashr.R
  #https://github.com/willwerscheid/flashier/blob/master/R/flash.R

  if(flashier_mix$n_factors == 0)
  {
    n = nrow(X)
    m = ncol(X)
    return(retDat(matrix(0,m,K), matrix(0,n,K), matrix(0,n,m),X, 0))
  }

  if(filter){
    factors <- filterFactors(flashier_mix, 0.01)
  }else
  {
    factors <- flashier_mix$loadings.pm[[2]]
  }

  loadings <- flashier_mix$loadings.pm[[1]]
  return(retDat(factors, loadings, fitted(flashier_mix),X, flashier_mix$pve))

}

#Helper functions
#filter
filterFactors <- function(flash_results, thresh = 0.01)
{
  keep <- flash_results$loadings.lfsr[[2]] > thresh
  filtered <- flash_results$loadings.pm[[2]]
  for(i in 1:ncol(filtered))
  {
    filtered[,i][keep[,i]] <- 0
  }
  return(filtered)
}

tune_cor_PMA <- function(X, rankK, su, sv){
  N = dim(X)[1]
  Tn = dim(X)[2]
  random_idx = arrayInd(sample(N*Tn,N*Tn, replace =F),dim(X))
  gap = floor(N / rankK)
  mres = 0
  for(i in seq(1,rankK)){
    rd_idx = random_idx[gap*(i-1)+1:gap*i, ]
    X_masked = X
    X_masked[rd_idx] = NA
    pmas = PMD(as.matrix(X_masked), K=rankK, sumabs = NULL, sumabsu = su, sumabsv = sv)
    X_hat = pmas$u %*% diag(pmas$d) %*% t(pmas$v)
    mres  = mres + sum((X_hat - X)[rd_idx]**2, na.rm = T) / (length(rd_idx) - sum(is.na((X_hat - X)[rd_idx])))
  }
  return(mres / rankK)
}

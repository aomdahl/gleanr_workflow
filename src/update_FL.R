################################################################################################################################
#Modified 4/26 by Ashton Omdahl
## Copied over from another function since this wasn't working originally.
################################################################################################################################

#helpful code for debugging
if(FALSE)
{
  #source("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/quickLoadData.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/quickLoadData.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_F.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_L.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/compute_obj.R")
  source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/plot_functions.R")
  all <- quickLoadFactorization("Z", "MARCC")
  X <- all$X
  W <- all$W
  option <- all$option
  option$traitSpecificVar <- TRUE
  option$parallel <- FALSE
  
  #subset for faster running....
  X <- X[1:1000, 1:10]
  W <- W[1:1000, 1:10]
  option$K <- 5
}

Update_FL <- function(X, W, option, preF = NULL, preL = NULL){
  # number of features - to avoid using T in R
  D = ncol(X)
  tStart0 = Sys.time()
  FactorM = NULL 
  #Random initialization
  if(!option[['preinitialize']])
  {
    FactorM   = matrix(runif(D*(option[['K']] - 1)), nrow = D);
    #Ashton added in- option to include column of 1s
    if(option[['ones_plain']])
    {
      message("Initializing ubiquitous factor with all ones...")
      FactorM = cbind(rep(1, D), FactorM);
    } else if(option[['ones_mixed_std']]) {
      message("Using a standard ubiq factor")
      cor_struct <- cor(X)
      #set up directions (-1 or 1) based on most common direction across the board.
      negatives <- apply(cor_struct, 1, function(x) sum((x < 0)))
      positives <- apply(cor_struct, 1, function(x) sum(x > 0) - 1) #-1 because every one is positively correlated with itself.
      best_ref <- positives >= negatives #this provides the most common reference
      #1 accept the most common reference as the directions- just assume this is the dominating effect
      ones <- ifelse(best_ref, 1, -1)
    } else if(option[['ones_mixed_ref']]) {
      cor_struct <- cor(X)
      #set up directions (-1 or 1) based on most common direction across the board.
      negatives <- apply(cor_struct, 1, function(x) sum((x < 0)))
      positives <- apply(cor_struct, 1, function(x) sum(x > 0) - 1) #-1 because every one is positively correlated with itself.
      best_ref <- positives >= negatives #this provides the most common reference
      #2 find some trait that matches this as closely as possible, let that be the reference and use the directions relative to that one.
      cor_ref <- cor_struct >= 0
      reference_trait <- which.max(apply(cor_ref, 2, function(x) sum(x == best_ref)))
      ones <- ifelse(cor_struct[,reference_trait] >= 0, 1, -1)
    }	else if(option[['ones_eigenvect']]) {
      message("1st column based on direction of svd of cor")
      cor_struct <- cor(X)
      svd <- svd(cor_struct, nu = 1) #fortunately its symmetric, so  U and V are the same here!
      ones <- sign(svd$u)
    } else {
      #FactorM   = matrix(runif(D*(option[['K']])), nrow = D);
      ones = matrix(runif(D*(option[['K']])), nrow = D)[,1]
    }
    FactorM = cbind(ones, FactorM);
  } else #you pre-provide the first F.
  {
    FactorM   = preF;
  }
  if(option$posF)
  {
    message("Performing semi-non-negative factorization today....")
    FactorM = abs(FactorM)
  }

  # First round of optimization
  message('Start optimization ...')
  message(paste0('K = ', (option[['K']]), '; alpha1 = ', (option[['alpha1']]),'; lambda1 = ', (option[['lambda1']])));
  
  #Need to adjust for first run- no reweighting (obvi.)
  og_option <- option[['reweighted']]
  option[['reweighted']] <- FALSE
  
  if(option[['parallel']]) #This is not working at all. Can't tell you why. But its not. Need to spend some time debugging at some point.
  {
    #print("Fitting L")	
    L = fit_L_parallel(X, W, FactorM, option, formerL = preL); #preL is by default Null, unless yo specify!
  } else
  {
    L = fit_L(X, W, FactorM, option, formerL = preL);
  }
  
  option[['reweighted']] <- og_option
  objective = c(NA, compute_obj(X, W, L, FactorM, option));
  objective_change = c(1, 1);
  old_change = 1;
  F_old = FactorM; 
  
  # if L is empty, stop
  non_empty_l = which(apply(L, 2, function(x) sum(x!=0)) > 0)
  if(length(non_empty_l) == 0){
    message('Finished');
    message('L is completely empty, alpha1 too large')
    FactorM = NULL;
    F_sparsity = 1;
    L_sparsity = 1;
    factor_corr = 1;
    Nfactor = 0;
    return()
  }
  
  trait.var <- matrix(NA, option[['iter']], ncol(X))

  for (iii in seq(1, option[['iter']])){
    message(paste0("Currently on iteration ", iii))
    ## update F
    #message("fitting F...")
    if(iii == 1)
    {
      og_option <- option[['reweighted']]
      option[['reweighted']] <- FALSE
      FactorM = fit_F(X, W, L, option, FactorM);
      option[['reweighted']] <- og_option
    } else {
      FactorM = fit_F(X, W, L, option, FactorM); #by iter 3 really slows down, due to the L1 requirements. Yea this won't do....
    }
    #get the factor specific variance....
    if(option$traitSpecificVar)
    {
      trait.var[iii,] <- FactorM$r.v
      FactorM <- FactorM$mat
    }
    # if number of factors decrease because of empty factor, the change in ||F||_F = 100
    if(ncol(FactorM) != ncol(F_old)){
      F_change = 100
    }else{
      F_change = norm(FactorM - F_old, 'F') / ncol(FactorM)
    }
    F_old = FactorM;
    non_empty_f = which(apply(FactorM, 2, function(x) sum(x!=0)) > 0)
    if(length(non_empty_f) == 0){
      message('Finished');
      message('F is completely empty, lambda1 too large')
      F_sparsity = 1;
      L_sparsity = 1;
      factor_corr = 1;
      Nfactor = 0;
      return();
    }
    colnames(FactorM) = seq(1, ncol(FactorM));
    
    ## update L
    #message("Fitting L....")
    if(option[['parallel']])
    {
      L  = fit_L_parallel(X, W, FactorM, option, L);
    }else{
      L  = fit_L(X, W, FactorM, option, L, r.v = trait.var[iii,]); #the l1 one requirement is making things tough here.....
    }

    # if L is empty, stop
    non_empty_l = which(apply(L, 2, function(x) sum(x!=0)) > 0)
    if(length(non_empty_l) == 0){
      message('Finished');
      message('L is completely empty, alpha1 too large')
      FactorM = NULL;
      F_sparsity = 1;
      L_sparsity = 1;
      factor_corr = 1;
      Nfactor = 0;
      return()
    }
    colnames(L) = seq(1, ncol(L));
    
    # align the two matrices	
    non_empty = intersect(non_empty_l,non_empty_f);
    L = as.matrix(as.data.frame(L[, non_empty]));
    FactorM  = as.matrix(as.data.frame(FactorM[, non_empty]));
    
    # collect sparsity in L and F
    L_sparsity = sum(abs(L) < 1e-5) / ncol(L) / nrow(L);
    F_sparsity = sum(abs(FactorM) < 1e-5) / ncol(FactorM) / nrow(FactorM);
    Nfactor = length(non_empty);
    
    # change in the objective -- should always be negative
    obj_updated = compute_obj(X, W, L, FactorM, option);
    obj_change = obj_updated - objective[length(objective)];
    objective = c(objective, obj_updated);
    objective_change = c(objective_change, obj_change);
    
    if(option[['disp']]){
      cat('\n')
      message(paste0('Iter', iii, ':'))
      message(paste0('Objective change = ', obj_change))
      message(paste0('Frobenius norm of (updated factor matrix - previous factor matrix) / number of factors  = ', F_change));
      message(paste0('Loading Sparsity = ', L_sparsity, '; Factor sparsity = ', F_sparsity, '; ', Nfactor, ' factors remain')); 
      cat('\n')
    }
    # converge if: 1). Change of the values in factor matrix is small. ie. The factor matrix is stable. 2). Change in the objective function becomes small; 3). reached maximum number of iterations
    if(option[['convF']] >= F_change){
     
      message(paste0('Factor matrix converges at iteration ', iii));
      break
    }
    
    oc1 = abs(objective_change[length(objective_change)])
    if(option[['convO']] >= oc1){
      message("Ojbective function change threshold achieved!")
      message(paste0('Objective function converges at iteration ', iii));
      break
    }
    if(iii == option[['iter']]){
      message('Reached maximum iteration.');
      break
    }
  }
  retlist <- list("F" = FactorM, "L" = L, "L_sparsity" = L_sparsity, "F_sparsity" = F_sparsity, "K" = Nfactor, "obj" = objective, "study_var" = trait.var)
  tEnd0 = Sys.time()
  cat('\n')
  message('Total time used for optimization: ');
  message(tEnd0 - tStart0);
  cat('\n')
  # return L, F, sparsity in L and F, number of factors -- could be different from K!
    return(retlist)
}

#9/07 testing notes
#tried regular, works
#trying now with factor-specific variance, shall see.
if(FALSE)
{

  #Get trait sample sizes:
  samp.size <- fread("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.n.tsv")
 ss <- data.frame(t(samp.size[1,2:11])) %>% mutate("names" = n[1:10])
 names(ss) <- c('count', "n")
 ss %>% arrange(count)
    #ran
    t2 <- Update_FL(X,W, option) #gets caught on iteration 5 here.... super weird.
    library(ggplot2)
    ndf <- data.frame("n" = names, "name" = paste0("X", 1:55) )
    res <- data.frame(t2[[7]]) %>% mutate("iter" = 1:30) %>% pivot_longer(cols = paste0("X", 1:10)) %>% merge(., ndf) %>%
      merge(., ss) %>% mutate("label" = paste0(n, "_", count))
    ggplot(res, aes(x = iter, y  = value, color = n)) + geom_line()
  
    
    ggplot(res, aes(x = iter, y  = value, color = label)) + geom_line(size = res$count/max(res$count))
    plotFactors(t2[[1]], trait_names = n[1:10], title = "mini")
    
    
    option$alpha1 <- 5
    option$lambda1 <- 5
    option$iter <- 10
    t3 <- Update_FL(X,W, option)
    res <- data.frame(t3[[7]]) %>% mutate("iter" = 1:10) %>% pivot_longer(cols = paste0("X", 1:10)) %>% merge(.,ndf) %>%
      merge(., ss) %>% mutate("label" = paste0(n, "_", count))
    ggplot(res, aes(x = iter, y  = value, color = label)) + geom_line(size = res$count/max(res$count))
    plotFactors(t3[[1]], trait_names = n[1:10], title = "mini")
    
    #Account for heritability 
    
    
    #Trying full,
    load("/work-zfs/abattle4/ashton/snp_networks/scratch/trait_standard_error_investigation/full.run.sept.RData")
    #X is full X
    option$alpha1 <- 20
    option$lambda1 <- 20
    option$iter <- 20
    samp.size <- fread("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.n.tsv")
    ss <- data.frame(t(samp.size[1,2:ncol(samp.size)])) %>% mutate("names" = n)
    names(ss) <- c('count', "n")
    ndf <- data.frame("n" = names, "name" = paste0("X", 1:55) )
    t.full <- Update_FL(X,W, option)
    res <- data.frame(t.full[[7]]) %>% mutate("iter" = 1:20) %>% pivot_longer(cols = paste0("X", 1:55)) %>% merge(.,ndf) %>%
      merge(., ss) %>% mutate("label" = paste0(n, "_", count))
    ggplot(res, aes(x = iter, y  = value, color = label)) + geom_line(size = res$count/max(res$count))  + theme(legend.position = "none")
    
    #I know this is so messy... diagnosing now.
  lowvar <- which(t.full[[7]][20,] < 0.01)
  names[lowvar] #I speculate these are the ones we really perform well on...
  plotFactors(t.full[[1]], trait_names = n, title = "all")
  

  #based on the factors, my predicted list
  
  #it looks like in general, they are mostly the top ones
  #look at reconstruction error
  X_hat <- t.full[[2]] %*% t(t.full[[1]])
  mse <- apply((X-X_hat)^2, 2, mean)
  plot(mse)
  which.min(mse)
  recon_order <-names[order(mse)]
  which(recon_order %in% names[lowvar])
  #exactly those 15.
  #So those that get low variance are those that we predict very well on.....
  #Let's look at the opposite end
  recon_order <-names[order(-mse)]
  highvar
  which(recon_order %in% names[lowvar])
  
  #plot of residual variance estimates bs actual
  e.var <- apply((X-X_hat),2,var)
  plot(apply((X-X_hat),2,var), t.full[[7]][20,], xlab = "Overall residual variance", ylab = "Output residual variance")
  abline(a = 0,b =1, col = "blue")
  diff <- abs(e.var -t.full[[7]][20,])
  mismatch <- which(abs(e.var -t.full[[7]][20,]) > 0.01)
  
  #Try it with a reduced regularization
  
  option$alpha1 <- 10
  option$lambda1 <- 10
  option$K <- 10
  option$iter <- 20
  samp.size <- fread("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.n.tsv")
 
  t.full.10 <- Update_FL(X,W, option)
  
}

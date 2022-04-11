pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr, penalized, cowplot, parallel, doParallel, flashr)
#source("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/sn-spMF/gwas_spMF/src/fit_F.R")
#source("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/sn-spMF/gwas_spMF/src/update_FL.R")
#source("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/sn-spMF/gwas_spMF/src/fit_L.R")
#source("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/sn-spMF/gwas_spMF/src/plot_functions.R")
#source('/Users/ashton/Documents/JHU/Research/LocalData/snp_network/sn-spMF/gwas_spMF/src/compute_obj.R')
quickLoadFactorization <- function(error_type, local)
{
  if(local == "MARCC")
  {
    dir = "/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/"
    names <- scan("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.names.tsv", what = character())
    n <- names
    samp.size <- fread("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.n.tsv")
  }else
  {
    dir = "~/Documents/JHU/Research/LocalData/snp_network/"
    names <- scan("~/Documents/JHU/Research/LocalData/snp_network/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1.names.tsv", what = character())
    n <- scan("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/matrix_sparsification/seed2_thresh0.9_h2-0.1.names.tsv", what = character())
  }
  z_scores <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.z.tsv")) %>% drop_na() %>% arrange(ids)
  pvals <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.p.tsv")) %>% drop_na() %>% arrange(ids)
  pvals <- as.matrix(pvals %>% select(-ids))
  z <- as.matrix(z_scores %>% select(-ids))
  any(duplicated(z_scores$ids))
  samp.sizes <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.z.tsv")) %>% drop_na() %>% arrange(ids)
  rsid_map <- fread(paste0(dir, "/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.pruned_rsids.txt"), header =FALSE) %>% rename("ids" = V1, "rsids" = V2) %>% 
    filter(ids %in% z_scores$ids) %>% .[!duplicated(.$ids),] %>% arrange(ids)
stopifnot(!any(rsid_map$ids != z_scores$ids))
  
  W_se <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.se.tsv")) %>% 
    drop_na() %>% filter(ids %in% z_scores$ids) %>% arrange(ids) %>% select(-ids) %>% as.matrix()
  W_se <- 1/ W_se
  
  W_maf <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.maf.tsv")) %>% 
    drop_na() %>% filter(ids %in% z_scores$ids) %>% arrange(ids) %>% select(-ids) %>% as.matrix()
  
  W_varmaf <- 1/matrix(apply(W_maf, 2, function(x) 2*x*(1-x)), nrow = nrow(W_maf), ncol = ncol(W_maf))
  
  betas <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.beta.tsv")) %>% drop_na() %>% filter(ids %in% z_scores$ids) %>% arrange(ids) %>% select(-ids) %>% as.matrix() 
  all_ids <- fread(paste0(dir, "/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.beta.tsv")) %>% drop_na() %>% filter(ids %in% z_scores$ids) %>% arrange(ids) %>% select(ids)
  W_ones <- matrix(1,nrow = nrow(W_maf), ncol = ncol(W_maf))
  #id like to have the sample sizes too....  
  if(error_type == "B_SE")
  {
    ret_dat <- betas
    ret_err <- W_se
  } else if(error_type == "B_MAF"){
    ret_dat <- betas
    ret_err <- W_varmaf
    
  }else{
    ret_dat <- z
    ret_err <- matrix(1, nrow(z), ncol(z))
  }
  option <- list()
  option[['K']] <- 15
  option[['alpha1']] <- 1
  option[['lambda1']] <- 1
  option[['iter']] <- 30
  option[['ones_mixed_std']] <- FALSE
  option[['ones_mixed_ref']] <- FALSE
  option[['ones']] <- FALSE
  option[["fixed_ubiq"]] <- FALSE
  option[['disp']] <- FALSE
  option[['convF']] <- 0
  option[['convO']] <- 0
  option[['ones_mixed_std']] <- FALSE
  option[["ones_mixed_ref"]] <- FALSE
  option[['ones_eigenvect']] <- TRUE
  option[['ones_plain']] <- FALSE
  option[['reweighted']] <- FALSE
  option[["glmnet"]] <- FALSE
  option[["parallel"]] <- TRUE
  option[["ridge_L"]] <- FALSE
  option[["preinitialize"]] <- FALSE
  option[['fastReg']] <- FALSE
  option$traitSpecificVar <- TRUE
  #Store stats here
  return(list("option" = option, "X" = ret_dat, "W" = ret_err, "names" = n, "pvals" = pvals, "vars" = rsid_map, "MAF"=W_maf))
}

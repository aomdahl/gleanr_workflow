#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pacman::p_load(magrittr, tidyr, dplyr, ggplot2, data.table)
BuildYAML <- function(row)
{
  ret.list <- list()
  dir = "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/"
  #factors
  ret.list[["factors"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$V, Category == "factors") %>% select(Directory, File)), collapse = "/"))
  #loadings
  ret.list[["loadings"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$U, Category == "loadings") %>% select(Directory, File)), collapse = "/"))
  #MAF
  ret.list[["maf"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$MAF, Category == "MAF") %>% select(Directory, File)), collapse = "/"))
  #K
  ret.list[["K"]] <- row$K
  #iter
  ret.list[["iter"]] <- row$NITER
  #sample size and overlap
  ret.list[["samp_overlap"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$N_o, Category == "samp_overlap") %>% select(Directory, File)), collapse = "/")) %>%
    gsub(x=., pattern = "_\\*_", replacement = paste0("_", as.character(row$N), "_"))

  #phenotype correlation
  ret.list[["pheno_corr"]] <- paste0(dir, paste0(unlist(filter(conv, Name == row$pheno_overlap, Category == "pheno_corr") %>% select(Directory, File)), collapse = "/"))

  #test_methods
  ret.list[["test_methods"]] <- row$test_methods
  #last things (fixed)
  ret.list[["noise_scaler"]] <- 1
  ret.list[["bic_param"]] <- "sklearn_eBIC"
  ret.list[["init"]] <- "V"
  if(is.na(row$HERIT))
  {
  ret.list[["herit_scaler"]] <- "continuous"

  }else
 {
	 ret.list[["herit_scaler"]] <- row$HERIT
  }
if(is.na(row$SHRINKAGE))
  {
	ret.list[["covar_shrinkage"]] <- -1
  }else
  {
	  ret.list[["covar_shrinkage"]] <- row$SHRINKAGE
  }


  data.frame("names" = names(ret.list), "vals" = unlist(ret.list))
}


#paramfiles <- "./OneDrive - Johns Hopkins/Research_Ashton_MacBook_Pro/snp_network/simulation_to_build.csv"
#conv.files <- "./OneDrive - Johns Hopkins/Research_Ashton_MacBook_Pro/snp_network/names_to_files_sims.csv"
#conv.files <- "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/setting_files/path_reference.csv"
conv.files <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/manuscript_analyses/simulations/2_factor_based/setting_files/path_reference.csv"
#paramfiles <- "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/setting_files/v1_u1_sims.csv"
#opath <- "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/25_round_sims_april23/"
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  message("Usage as follows:")
  message("parameter_list_file.csv output_path <optional> -c")
  message("-c : this will print out the commands to run the simulation at the specified settings at the destination path specified with -c")
  stop("At least two arguments must be supplied (input file path and output directory).n", call.=FALSE)

}
paramfiles <- args[1]
opath <- args[2]
param.setting <- fread(paramfiles)


param.setting <- param.setting %>% mutate("fnames" = paste0(V,"_", U, "_MAF-", MAF, "_N-",N, "_RHO-", pheno_overlap, "_No-", N_o, ".yml"))


conv <-fread(conv.files)
cats <- c("factors", "loadings", "maf", "K", "iter", "samp_overlap","pheno_corr", "test_methods","noise_scaler","bic_param", "init", "shrinkage")
if (length(args)==4) {
cat("mkdir -p ", args[4], "\n")
  }
for(i in 1:nrow(param.setting))
{

  if (length(args)==4) {
    if(args[3] == "-c")
    {
      cat("bash src/runSimulation.sh ", paste0(opath, '/',param.setting[i,]$fnames), paste0(args[4], gsub(x=param.setting[i,]$fnames, replacement = "", pattern = ".yml"),"/"), "\n")
    }
  }else
  {
    message("Currently running ", param.setting[i,]$fnames)
  }

  tab <- BuildYAML(param.setting[i,])
  write.table(tab, file = paste0(opath, "/", param.setting[i,]$fnames), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
}



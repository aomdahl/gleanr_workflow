############################################################################################################
# Ashton Omdahl, November 2024
### Miscellaneous supplemental figures:
############################################################################################################
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")
source("/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/src/get_pLI.R")

pacman::p_load(magrittr, dplyr, ggplot2, data.table, ggsignif,RColorBrewer)
# Supplementary table 1- phenotypes included in UKBB/Finngen analysis:
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")

#Update this to something reflect
png('/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/suppl_fig_gleanr_convergence.png')
plot(sapply(ret$each.matrix.obj, function(x) x$To), pch=19, col="skyblue", xlab="ALS Iteration",ylab="Objective function")
dev.off()

#For figure of enrichments, see scratch/factor_interpretatin/analysis/blood_factor_enrichr_results.Rmd``


#relationship between K init and K
path <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/"
data <- list.files(path="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/", pattern ="_bic_dat.RData" )
curr <- NULL

for(d in data)
{
  min=0
  load(paste0(path, d))
  k <- gsub(x=d, pattern = "final_gleaner_attempt_K", replacement = "") %>% gsub(x=., pattern ="_bic_dat.RData", replacement = "" )
  if(any(apply(bic.dat$optimal.v,2,function(x) all(x==0))))
  {
    min <- sum(apply(bic.dat$optimal.v,2,function(x) all(x==0)))
  }
  curr <- rbind(curr,c(k,ncol(bic.dat$optimal.v)-min))
}

curr <- data.frame(curr) %>% set_colnames(c("K_init", "K_final"))
plot(curr$K_init, curr$K_final)
#not great. Maybe drop it...

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/conservative_1e-5/covar_influence/BIC-sklearn_eBIC_K-GRID_K_search.RData")
grid.search.record$curr_runs
ncol(grid.search.record$test_dat$`25%`$optimal.v)
plot(grid.search.record$curr_runs$Kinit, grid.search.record$curr_runs$K_end)

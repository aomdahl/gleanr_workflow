getIterDetails <- function(k_search_obj, source)
{
  data.frame("N"=nrow(k_search_obj$test_dat$`25%`$optimal.u), "M"=nrow(k_search_obj$test_dat$`25%`$optimal.v),"K_init"= k_search_obj$curr_runs$Kinit,"n_updates_l"=sapply(k_search_obj$test_dat, function(x) length(x$rec.dat$bic.l)),
             "n_updates_all"=sapply(k_search_obj$test_dat, function(x) length(x$rec.dat$bic_sum)),
             "source"=source)
}

#Global run, big one
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K68_bic_dat.RData")
stat.mat <- c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=68, "n_updates_l"=length(bic.dat$rec.dat$bic.l), "n_updates_all"=length(bic.dat$rec.dat$bic_sum), "source"="real_data")

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K102_bic_dat.RData")
stat.mat <- rbind(stat.mat, c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=102, "n_updates_l"=length(bic.dat$rec.dat$bic.l), n_updates_all=length(bic.dat$rec.dat$bic_sum), "source"="real_data"))

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K119_bic_dat.RData")
stat.mat <- rbind(stat.mat, c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=119, "n_updates_l"=length(bic.dat$rec.dat$bic.l), n_updates_all=length(bic.dat$rec.dat$bic_sum), "source"="real_data"))

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K34_bic_dat.RData")
stat.mat <- rbind(stat.mat, c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=34, "n_updates_l"=length(bic.dat$rec.dat$bic.l), n_updates_all=length(bic.dat$rec.dat$bic_sum), "source"="real_data"))

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K127_bic_dat.RData")
stat.mat <- rbind(stat.mat, c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=34, "n_updates_l"=length(bic.dat$rec.dat$bic.l), n_updates_all=length(bic.dat$rec.dat$bic_sum), "source"="real_data"))

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K136_bic_dat.RData")
stat.mat <- rbind(stat.mat, c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=34, "n_updates_l"=length(bic.dat$rec.dat$bic.l), n_updates_all=length(bic.dat$rec.dat$bic_sum), "source"="real_data"))

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K134_bic_dat.RData")
stat.mat <- rbind(stat.mat, c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=34, "n_updates_l"=length(bic.dat$rec.dat$bic.l), n_updates_all=length(bic.dat$rec.dat$bic_sum), "source"="real_data"))



load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K130_bic_dat.RData")
stat.mat <- rbind(stat.mat, c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=34, "n_updates_l"=length(bic.dat$rec.dat$bic.l), n_updates_all=length(bic.dat$rec.dat$bic_sum), "source"="real_data"))


#3Now on the finngen/ukbb set

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/covar_influence/BIC-sklearn_eBIC_K-41_bic_dat.RData")
stat.mat <- rbind(stat.mat, c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=41, "n_updates_l"=length(bic.dat$rec.dat$bic.l), n_updates_all=length(bic.dat$rec.dat$bic_sum), "source"="real_data"))

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/covar_influence/BIC-sklearn_eBIC_K-GRID_K_search.RData")
stat.mat <- rbind(stat.mat, getIterDetails(grid.search.record, "real_data"))
  

##UKBB ste

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2//conservative_1e-5/covar_influence/BIC-sklearn_eBIC_K-41_bic_dat.RData")
stat.mat <- rbind(stat.mat, c("N"=nrow(bic.dat$optimal.u), "M"=nrow(bic.dat$optimal.v), "K_init"=41, "n_updates_l"=length(bic.dat$rec.dat$bic.l), n_updates_all=length(bic.dat$rec.dat$bic_sum), "source"="real_data"))

load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/conservative_1e-5/covar_influence/BIC-sklearn_eBIC_K-GRID_0_K_search.RData")
stat.mat <- rbind(stat.mat,getIterDetails(grid.search.record, "real_data"))

#Some from the sims
simdir <- "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/"
load(paste0(simdir, '/2b_overlap/V103_U103_MAF-mix_eur_N-50000_RHO-2b_mid_mixed_p_No-2b_high_no/factorization_results/sim10.GLEANER_glmnet_K_search.RData'))
stat.mat <- rbind(stat.mat,getIterDetails(grid.search.record, "sims"))


load(paste0(simdir, '/1b_overlap/V101_U103_MAF-mix_eur_N-10000_RHO-1b_high_mixed_p_No-1b_high_no/factorization_results/sim5.GLEANER_glmnet_K_search.RData'))
stat.mat <- rbind(stat.mat,getIterDetails(grid.search.record, "sims"))


load(paste0(simdir, '/no_overlap/V102_U101_MAF-mix_eur_N-2e+05_RHO-none_No-none/factorization_results/sim2.GLEANER_glmnet_K_search.RData'))
stat.mat <- rbind(stat.mat,getIterDetails(grid.search.record, "sims"))



library(ggplot2)
library(dplyr)
library(magrittr)
stat.mat <- data.frame(stat.mat) %>% mutate_at(all_of(c("N", "M", "n_updates_l","n_updates_all", "K_init")), as.numeric)
#Not a very strong trend- appeasr to decline with Kinit, M and N
ggplot(stat.mat, aes(x=as.factor(N),y=n_updates_l, color=K_init, size=M)) + geom_jitter() + theme_bw() + facet_wrap(~source, scales="free_x") + 
  xlab("# SNPs (N)")+ ylab("# alternating updates")
ggsave(file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/supp_10_bic.png",height=5, width=10)




suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("magrittr"))
suppressMessages(library("dplyr"))
suppressMessages(library(stringr))
suppressMessages(library("combinat"))
source("../../src/matrix_similarity_scoring/evaluateSimR2.R")
#source("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/src/evaluateSimR2.R")
#source("/scratch16/abattle4/ashton/snp_networks/scratch/finngen_v_ukbb/src/matrix_comparison_utils.R")
source("../../src/matrix_similarity_scoring/matrix_comparison_utils.R")
library("gleanr")
#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/R/sparsity_scaler.R")

suppressMessages(library(NatParksPalettes))
option_list <- list(
  make_option(c("-o", "--output"), default = '', type = "character",
              help="Output dir + handle"),
  make_option(c("-s", "--sim_path"), default = '', type = "character",
              help="File containing all of the paths of simulations you wish to go into this plot."),
  make_option(c("-p", "--plot"), default = FALSE, type = "logical", action = "store_true",
              help="specify this if you want to plot."),
  make_option(c("-y", "--yaml"), default = '', type = "character",
              help="Path to the settings file."),
  make_option(c("-w", "--whitened"), default = FALSE, type = "logical",
              help="Include whitened", action = "store_true"),
  make_option(c("-n", "--scale_data"), default = FALSE, type = "logical",
              help="Scale true and predicted V/U to have unit norm.", action = "store_true"),
  make_option(c("--sim_raw"), default = "", type = "character",
              help="Path to the actual simulation files, not the tests")
  
)

#debug interactive

t = c("--output=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.3_n100gwasMF_optim_covar_high_fixed_1/factorization_results/summary",
  "--yaml=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/V1_U1_maf0.3_n100.high_covar_fixed.yml",
  "--sim_path=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.3_n100gwasMF_optim_covar_high_fixed_1/factorization_results/",  "--scale_data")
#V1_U1_maf0.01_n100.yml.gwasMF_only.yml
t = c("--output=//scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V7_U7_mafmixed_n50000.no_covar_cont_scaling//factorization_results/summary",
      "--yaml=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/V7_U7_mafmixed_n50000.no_covar_cont_scaling.yml",
      "--sim_path=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V7_U7_mafmixed_n50000.no_covar_cont_scaling//factorization_results/",  "--scale_data")

finalpath="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/1b_overlap/"
t = c(paste0("--output=", finalpath, "V101_U102_MAF-mix_eur_N-50000_RHO-1b_high_mixed_p_No-1b_high_no/factorization_results/summary"),
      "--yaml=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/final_sims_june_2024/1b_overlap/V101_U102_MAF-mix_eur_N-50000_RHO-1b_high_mixed_p_No-1b_high_no.yml",
      "--scale_data", paste0("--sim_path=", finalpath, "V101_U102_MAF-mix_eur_N-50000_RHO-1b_high_mixed_p_No-1b_high_no/factorization_results/"))

#args <- parse_args(OptionParser(option_list=option_list), args = t)

args <- parse_args(OptionParser(option_list=option_list))#, args = t)

#                   args = c("--input=/Users/aomdahl/Library/CloudStorage/OneDrive-JohnsHopkins/Research_Ashton_MacBook_Pro/snp_network/scratch/cohort_overlap_exploration/babytest1"))
if(args$sim_raw =="")
{
  message("setting sim raw directory as one up from sim_path")
  args$sim_raw = gsub(pattern="factorization_results\\/",x=args$sim_path, replacement="")
}
#Rscript src/evaluateSims.R --output simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/factorization_results/ -p -y  simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/udler_1_simple_no-corr_no-overlap.yml-s simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/factorization_results/
#Organized by directory per sim.
#corresponding script is in ldsc_all_traits/src/factMetaAssessment
#handles file, i.e. pma,pmd,flashr....
#$yaml
#[1] "./simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/udler_1_simple_no-corr_no-overlap.yml"
#args$sim.path <- "./simulating_factors/udler_based_500/udler_1_simple_no-corr_no-overlap/factorization_results/"
#get the names of the methods you want

yml <- read.table(args$yaml, sep = ",")
methods.run <- ((yml %>% filter(V1 == "test_methods"))$V2 %>% strsplit(.,":" ))[[1]]
true.loadings <- as.matrix(fread((yml %>% filter(V1 == "loadings"))$V2))
true.factors <-  as.matrix(fread((yml %>% filter(V1 == "factors"))$V2))
true.betas <- as.matrix(fread(paste0(args$sim_raw, "/noise-free_effect_sizes.txt"))[,-1])
true.se <- as.matrix(fread(paste0(args$sim_raw, "/sim1.std_error.txt"))[,-1])
if(args$scale_data)
{
  true.loadings <- unitNorms(true.loadings)
  true.factors <- unitNorms(true.factors)
}
niter <- as.numeric((yml %>% filter(V1 == "iter"))$V2)
s <- args$sim_path
true.k <- ncol(true.factors)
#This will evaluate the simulations you desire. The assumption is that each individual simulation has its own directory, with a
#sim_loadings and sim_factors file, and then files corresponding to the predictved outcome for several different methods.
sim.performance <- NULL
f_i =1
nmethods = length(methods.run)
ref.list <- list()
pred.list <- list()
n_features=12
for(m in methods.run)
{
  message("Assessing simulations from ", m)
  r_performance <- matrix(NA, nrow = niter, ncol = n_features)
  for(i in 1:niter){
    lookup_id=paste0(m,"_", i)
    pred.list[[lookup_id]] <- list()
    ref.list[[lookup_id]] <- list("U"=true.loadings, "V"=true.factors)
    #print(i)
    #print(m)
    if(!file.exists(paste0(s, "/sim",i, ".", m, ".loadings.txt")) & i > 1)
      {
         message("WARNING: missing file: ", paste0(s, "/sim",i, ".", m, ".loadings.txt"))
         next
    }
    if(!file.exists(paste0(s, "/sim",i, ".", m, ".loadings.txt")) & i == 1)
       {
         message("WARNING: missing file: ",paste0(s, "/sim",i, ".", m, ".loadings.txt"))
         break;
    }
    #print(i)
    #print(paste0(s, "/sim",i, ".", m, ".loadings.txt"))
    pred.list[[lookup_id]][["U"]] <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".loadings.txt")))
    pred.list[[lookup_id]][["V"]] <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".factors.txt")))
    pred.list[[lookup_id]][["X_hat"]] <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".X-hat.txt")))
    if(args$scale_data)
    {
     # message("Scaling both true and loaded data for convenient comparison")
      if(!all(pred.list[[lookup_id]][["U"]] == 0))
      {
        pred.list[[lookup_id]][["U"]] <-unitNorms(pred.list[[lookup_id]][["U"]])
        pred.list[[lookup_id]][["V"]] <- unitNorms(pred.list[[lookup_id]][["V"]])
      }else
      {
        message("Unable to scale data; all elements 0'd out.")
      }

    }
    reconstruction <- fread(paste0(s, "/sim",i, ".", m, ".recon_error.txt"))

    #Get K
    ks <- ncol(pred.list[[lookup_id]][["V"]])
    sparsity.v <- matrixSparsity(pred.list[[lookup_id]][["V"]],initK = true.k, wrt.init = TRUE)
    sparsity.u <- matrixSparsity(pred.list[[lookup_id]][["U"]],initK = true.k, wrt.init = TRUE)
    if(all(pred.list[[lookup_id]][["U"]] == 0))
    {
      ks <- 0
    }

  #maybe have the option to scale be embedded in here?
    #The columns here are: c("method", 'R2_L','R2_F',"RSE", "R2_X", "K_out","sparsityV", "sparsityU", "iter")
    #compareModelMatricesComprehensive <- function(A,B, corr.type = "pearson", full.procrust=TRUE)
    #Test.findings
    #Old version
    yuan.style <- evaluate_error(true.loadings, true.factors, pred.list[[lookup_id]][["U"]], pred.list[[lookup_id]][["V"]])
    xhat_calc <- xhatFit(m, pred.list[[lookup_id]][["X_hat"]], true.betas, true.se) #return(list("rrmse"=rrmse(x_hat_scaled, x_true), "cor"=stackAndAssessCorr(x_hat_scaled, x_true)))
    r_performance[i,] <- c(evaluteFactorConstruction(true.loadings, true.factors, pred.list[[lookup_id]][["U"]], pred.list[[lookup_id]][["V"]],unit.scale = FALSE,null.return=TRUE),
                           reconstruction$Frobenius_norm[1], reconstruction$Correlation[1], yuan.style$U_r2, yuan.style$V_r2,"X_RRMSE"=xhat_calc$rrmse,"X_cor"=xhat_calc$cor,  ks,sparsity.v,sparsity.u, i)
    r_performance[i,n_features] <- i
    f_i = f_i+1
  }
  sim.performance <- rbind(sim.performance, data.frame("method" = m, r_performance))
}
##now the alternative version:
clean.dat <- list(makeTableOfComparisons(ref.list, pred.list, by_order = FALSE, filter_term = "*"))
updated.calcs <- extractAndTabularizeRunGeneric(clean.dat)

#Now, read in whitened if they are there....
if(args$whitened & FALSE)
{
  nmethods = length(methods.run) *2
  for(m in methods.run)
  {
    r_performance <- matrix(NA, nrow = niter, ncol = 5)
    for(i in 1:niter){
      if(!file.exists(paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt")) & i > 1)
      {
       message("WARNING: missing file: ", paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt"))
        next
      }
      if(!file.exists(paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt")) & i == 1)
         {
           message("WARNING: missing file: ", paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt"))
           break
      }
      pred.list[[lookup_id]][["U"]] <- fread(paste0(s, "/sim",i, ".", "whitened.", m, ".loadings.txt"))
      pred.list[[lookup_id]][["V"]] <- fread(paste0(s, "/sim",i, ".","whitened.", m, ".factors.txt"))
      if(args$scale_data)
      {
        message("Scaling both true and loaded data for convenient comparison")
        pred.loadings <- apply(pred.loadings, 2, function(x) x/norm(x, "2"))
        pred.list[[lookup_id]][["V"]] <- apply(pred.list[[lookup_id]][["V"]], 2, function(x) x/norm(x, "2"))
      }

      reconstruction <- fread(paste0(s, "/sim",i, ".","whitened.", m, ".recon_error.txt"))
      yuan.style <- evaluate_error(true.loadings, true.factors, pred.loadings, pred.factors)
      r_performance[i,] <- c(evaluteFactorConstruction(true.loadings, true.factors, pred.loadings, pred.list[[lookup_id]][["V"]]),
                             reconstruction$Frobenius_norm[1], reconstruction$Correlation[1], i)
      f_i = f_i+1
    }
    sim.performance <- rbind(sim.performance, data.frame("method" = paste0(m, "_whitened"), r_performance))
  }
}
colnames(sim.performance) <- c("method", 'R2_L','R2_F',"RSE", "R2_X","yuan_U", "yuan_V","X_RRMSE","X_cor", "K_out","sparsityV", "sparsityU", "iter")

#Combine the old and the new
#Updated August 27 to reflect the most recent match-based version
#full.sim.performance <- cbind(sim.performance, updated.calcs$V %>% dplyr::select(BIC, PRC, GlobalKappa,correlation, Procrustes_pearson, Euclidian_dist) %>% 
#                                set_colnames(c("runID_V", "PRC_V", "Kappa_V", "Corr_mine_V", "Procrust_pearson_V", "Distance_V")),
#updated.calcs$U %>% dplyr::select(BIC, PRC, GlobalKappa,correlation, Procrustes_pearson, Euclidian_dist) %>% set_colnames(c("runID_U", "PRC_U", "Kappa_U",  "Corr_mine_U", "Procrust_pearson_U", "Distance_U")),
#updated.calcs$Xhat_norms) %>% rename("Xhat_dist"=dist)

stopifnot(all(paste0(sim.performance$method, "_", sim.performance$iter) == updated.calcs$joint$BIC))
full.sim.performance <- cbind(sim.performance, updated.calcs$joint %>% 
                                dplyr::select(BIC, joint_cb_auc_v, joint_cb_kappa_v,joint_cb_pearson_v, joint_cb_kendall_v, joint_cb_rrmse_v,
                                              joint_pr_auc_v, joint_pr_kappa_v,joint_pr_pearson_v, joint_pr_kendall_v, joint_pr_rrmse_v,
                                              joint_cb_auc_u, joint_cb_kappa_u,joint_cb_pearson_u, joint_cb_kendall_u, joint_cb_rrmse_u,
                                              joint_pr_auc_u, joint_pr_kappa_u,joint_pr_pearson_u, joint_pr_kendall_u, joint_pr_rrmse_u),
                              updated.calcs$Xhat_norms %>% rename("Xhat_dist"=dist))

stopifnot(all(full.sim.performance$BIC == paste0(sim.performance$method, "_", sim.performance$iter)))
#if(all(full.sim.performance$runID_U == full.sim.performance$runID_V))
#{
#  if(all(paste0(full.sim.performance$method, "_", full.sim.performance$iter) == full.sim.performance$runID_V))
#  {
#    full.sim.performance <- full.sim.performance %>% select(-runID_V, -runID_U)
#  } else
#  {
#    message("ERROR: simulations don't line up in assessment")
#    quit()
#  }
#}

#if(any(full.sim.performance$R2_F > full.sim.performance$Procrust_pearson_V^2))
#{
#  message("WARNING: our R2 on V exceeds the procrustes version.")
#  print(full.sim.performance %>% filter(R2_F > Procrust_pearson_V^2))
#}
#if(any(full.sim.performance$R2_L > full.sim.performance$Procrust_pearson_U^2))
#{
#  message("WARNING: our R2 on V exceeds the procrustes version.")
#  print(full.sim.performance %>% filter(R2_L > Procrust_pearson_U^2))
#}

#Now plot it, if desired
write.table(full.sim.performance, file = paste0(args$output, ".tabular_performance.tsv"), quote = FALSE, row.names = FALSE)
if(args$plot & FALSE)
{
  #Plots of F and L
  suppressMessages(library(ggplot2))
  l.perf <- ggplot(sim.performance, aes(x = method, y = R2_L, fill = method)) + geom_boxplot() + ylab(expression(R^2)) +   ggtitle("Loadings") + theme_classic(15) +
    scale_fill_manual(values=natparks.pals("Arches", nmethods)) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank())
  f.perf <- ggplot(sim.performance, aes(x = method, y = R2_F, fill = method)) + geom_boxplot() + ylab(expression(R^2)) +
    ggtitle("Factors") + theme_classic(15) + scale_fill_manual(values=natparks.pals("Arches", nmethods))
  suppressMessages(library(cowplot))
  ggsave(plot_grid(plotlist = list(l.perf, f.perf), nrow = 2, ncol =1), filename = paste0(args$output, ".f_l_boxplot.png"))

    #plots of frobenious norm and overall R2
    frob <- ggplot(sim.performance, aes(x = method, y = RSE, fill = method)) + geom_boxplot() + ylab("Root squared error") +   ggtitle("Frobenius norm") + theme_classic(15) +
    scale_fill_manual(values=natparks.pals("Arches", (nmethods))) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank())
  r2x <- ggplot(sim.performance, aes(x = method, y = R2_X, fill = method)) + geom_boxplot() + ylab(expression(R^2)) +
    ggtitle("Overall correlation") + theme_classic(15) + scale_fill_manual(values=natparks.pals("Arches", nmethods))
  ggsave(plot_grid(plotlist = list(frob, r2x), nrow = 2, ncol =1), filename = paste0(args$output, ".recon_boxplot.png"))

  #K-plot
  ks <- ggplot(sim.performance, aes(x = method, y = K_out, fill = method)) + geom_boxplot() + ylab("Final K") +
    ggtitle("Predicted latent factors") + theme_classic(15) + scale_fill_manual(values=natparks.pals("Arches", nmethods))
  ggsave(ks, filename = paste0(args$output, ".ks_boxplot.png"))


}





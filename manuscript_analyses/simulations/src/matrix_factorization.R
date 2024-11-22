#!/usr/bin/Rscript --vanilla
pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, cowplot, flashr, ashr, optparse)
#Some useful functions
#From stack overflow
if("flashr" %in% rownames(installed.packages()) == FALSE) {
library(devtools)
install_github("stephenslab/flashr",build_vignettes = FALSE)
}


#Helper output funtions
plotMatrices <- function(dat, desc, trait_names, snp_ids){
	ggsave(plotFactors(dat$V,trait_names, paste0(desc, " Factors")), filename =paste0(args$outdir, desc, ".factors.png"), width = 11, height = 8)
  ggsave(plotLoadings(dat$U, snp_ids, paste0(desc, " Loadings ")), filename =paste0(args$outdir, desc, ".loadings.png"), width = 7, height = 10 )
}

writeMatrices <- function(dat, desc){
  write_tsv(data.frame(dat$V), paste0(args$outdir, desc,".factors.txt") )
  write_tsv(data.frame(dat$U), paste0(args$outdir, desc,".loadings.txt") )
  write_tsv(data.frame(dat$X_hat), paste0(args$outdir, desc,".X-hat.txt") )
}

#this will write out the
writeReconError <- function(dat, desc, X)
{
  frob <- norm(X - dat$X_hat, "F")
  #TODO- add a thing here where if everything is 0 or X is null, we just give an r of 0.
  if(all(as.vector(dat$X_hat) == 0))
  {
    r <- 0
  }else
  {
    r <- stats::cor(as.vector(X), as.vector(dat$X_hat))
  }
  
  write.table(x=data.frame("Frobenius_norm" = frob, "Correlation" = r),file=paste0(args$outdir, desc,".recon_error.txt"),row.names = FALSE, quote = FALSE )
}

#Also accounts for sparsity, cuz we can
writePVE <- function(dat, desc)
{
  if(length(dat$PVE) == 0)
  {
    w <- data.frame()
    write.table(x=w,file=paste0(args$outdir, desc,".pve.txt"),row.names = FALSE, quote = FALSE)
    return()
  }
  if(!is.null(dat$pve) & !is.na(dat$pve))
  {
    thresh = 1e-4
    if( (ncol(dat$V) != length(dat$PVE))  )
    {
      message("STOP HERE")

    }else    {
      w <- data.frame("Factor" = paste0("F", 1:length(dat$pve)), "PVE" = dat$pve, "Sparsity" = colSums(abs(dat$V) < thresh)/nrow(dat$V))
      write.table(x=w,file=paste0(args$outdir, desc,".pve.txt"),row.names = FALSE, quote = FALSE)
      return()
    }
  }
}

#writeAllOutputs(run.data, meth, run.data$X)
writeAllOutputs <- function(results, description, x_true)
{
  writeMatrices(results, description)
  writeReconError(results, description, x_true)
  writePVE(results, description)
}


#source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("../../src/plot_functions.R")
source("src/factorization_methods.R")
#source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/factorization_methods.R")
option_list <- list(
	make_option("--gwas_data", type = 'character', help = "path to gwas z-scores data; table of n snp rows x m study columns.", default = ""),
	make_option("--se_data", type = 'character', help = "path to gwas snp SEs data; table of n snp rows x m study columns."),
	make_option("--beta_data", type = 'character', help = "path to gwas betas; table of n snp rows x m study columns."),
	make_option("--trait_names", type = 'character', help = "List of trait names, given in same order as input matrix.", default = ""),
	make_option("--seed", type= "integer", default = 1, help = "Specify the seed for consistency."),
	make_option("--outdir", type  = "character", help = "Specify an output directory"),
	make_option('--test_run',type='logical',action='store_true',help='Do a test run of the program.', default = F),
	make_option('--all_run',type='logical',action='store_true',help='run all the factorization methods..', default = F),
	make_option('--only_run',type='character',help = "Specify only particular methods to run if you don't want all of them; separate by a comma.", default = ""),
	make_option('--single_method',type='character',help='Run just a single method, specify (i.e. ssvd, pca, flash, etc.)', default = ""),
	make_option('--no_plots',type='logical',help='Specify to not make any plots', default = FALSE, action = "store_true"),
	make_option('--K',type='integer',help='Specify a k to run', default = 10),
	make_option('--C',type='character',help='Specify a true noise covariance file to correct, if you want.', default = ""),
	make_option('--bic_var',type='character',help='Specify bic var type', default = "mle"),
	make_option('--step_scaling',type='logical',help='Specify bic var type', default = FALSE, action = "store_true"),
	make_option('--init_mat',type='character',help='Specify which matrix gets initialized', default = "V"),
	make_option('--WLgamma',type='numeric',help='SPecify the degree of WL shrinkage, if any.', default = -1),
	make_option('--no_column_id',type='logical',help='SPecify if no column id, default is there is.', default = TRUE),
	make_option('--z_scores',type='character',help='Path to z scores; n snps by m study columsn, a tsv file'),
	make_option(c("-n", "--n_samples"),type='character',help='Path to sample size per study file'),
	make_option('--debug',type='logical',action='store_true',help='Use if you are running on debug.', default = FALSE)
)
message("Beware the double negative argument no column id.")

#UKBB style sim
dir="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/ukbb_based_sims/seed1/se_scale_0.001/"
t=c(paste0("--se_data=", dir, "ukbb.sim1.se.tsv"),
    paste0("--beta_data=", dir, "ukbb.sim1.beta.tsv"),
    "--seed=1",
    paste0("--outdir=", dir, "sim1.manual"),
    "--only_run=PCA,PCA_chooseK,FLASH_SE,FLASH_SE_kroneker,GLEANER_glmnet,GLEANER_glmnet_noCovar", "--K=8", "--bic_var=mle",
    "--no_plots", "--init_mat=V")
#yuan Sim
t=c("--se_data=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/yuan_simulations/Input_tau100_seed1_W.txt",
    "--beta_data=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/yuan_simulations/Input_tau100_seed1_X.txt",
    "--seed=1",
    "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/yuan_simulations//simMANUAL.",
    "--only_run=GLEANER,PCA", "--K=5", "--bic_var=dev",
    "--no_plots", "--init_mat=V")
#
#####SANN_brent
t=c("--se_data=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.4_n50000.new_seed.SANN-Brent//sim3.std_error.txt",
    "--beta_data=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.4_n50000.new_seed.SANN-Brent//sim3.effect_sizes.txt",
    "--seed=100001", "--outdir=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.4_n50000.new_seed.SANN-Brent/factorization_results/sim3MaNUAL",
    "--only_run=GLEANER_glmnet", "--K=5", "--no_plots", "--bic_var=mle", "--init_mat=V",
    "--C=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.4_n50000.new_seed.SANN-Brent//sim3.c_matrix.txt")

t=c("--se_data=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/broad25Sims/V1_U1_MAF-mix_N-2500_RHO-1b_high_p_No-1b_high_no//sim10.std_error.txt",
    "--beta_data=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs//broad25Sims/V1_U1_MAF-mix_N-2500_RHO-1b_high_p_No-1b_high_no//sim10.effect_sizes.txt",
    "--seed=10", "--outdir=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/broad25Sims/V1_U1_MAF-mix_N-2500_RHO-1b_high_p_No-1b_high_no//factorization_results/sim10MaNUAL",
    "--only_run=GLEANER_noCovar", "--K=5", "--no_plots", "--bic_var=mle", "--init_mat=V",
    "--C=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/broad25Sims/V1_U1_MAF-mix_N-2500_RHO-1b_high_p_No-1b_high_no//sim10.c_matrix.txt")

#BFGS DOES WELL IN THIS CONTEXT
t=c("--se_data=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.4_n5000.high_covar_1block_SCALED_Wc_t//sim10.std_error.txt",
"--beta_data=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.4_n5000.high_covar_1block_SCALED_Wc_t//sim10.effect_sizes.txt",
"--seed=1", "--outdir=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.4_n5000.high_covar_1block_SCALED_Wc_t//factorization_results/sim10MaNUAL",
"--only_run=gwasMF_BIC", "--K=5", "--no_plots", "--bic_var=mle", "--init_mat=V",
"--C=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V1_U1_maf0.4_n5000.high_covar_1block_SCALED_Wc_t//sim10.c_matrix.txt")


#generic sim
t=c("--se_data=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V7_U7_mafmixed_n5000.COVAR_cont_scaling//sim8.std_error.txt",
    "--beta_data=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V7_U7_mafmixed_n5000.COVAR_cont_scaling/sim8.effect_sizes.txt",
    "--seed=1", "--outdir=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V7_U7_mafmixed_n5000.COVAR_cont_scaling/factorization_results/sim8MaNUAL",
    "--only_run=gwasMF_BIC", "--K=0", "--no_plots", "--bic_var=sklearn", "--init_mat=V",
    "--C=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs/V7_U7_mafmixed_n5000.COVAR_cont_scaling/sim8.c_matrix.txt")


## GWAS-BASED SIM, APRIL 2024- no overlap
simp="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/V103_U101_MAF-mix_eur_N-10000_RHO-none_No-none/"
t=c(paste0("--se_data=",simp,"seed4.replicate23.0.SE.csv"),
    paste0("--beta_data=",simp,"seed4.replicate23.0.BETA.csv"),
    "--seed=1",
    "--outdir=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/V103_U101_MAF-mix_eur_N-10000_RHO-none_No-none/",
    "--only_run=GLEANER_glmnet_noCovar", "--K=8",
    "--no_plots")

  # High overlap
  t=c(paste0("--se_data=",simp,"sim3.std_error.txt"),
      paste0("--beta_data=",simp,"sim3.effect_sizes.txt"),
      "--seed=1",
      "--outdir=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/sandbox/gleaner_test",
      "--only_run=GLEANER_glmnet,GLEANER_glmnet_noCovar", "--K=5",
      "--no_plots")
  
#Testing factorGo yikes

  
  over="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/1b_overlap_dev_grid/V101_U101_MAF-mix_eur_N-10000_RHO-1b_high_mixed_p_No-1b_high_no/"
  t=c(paste0("--se_data=",over,"sim1.std_error.txt"),
      paste0("--beta_data=",over,"sim1.effect_sizes.txt"),
      "--seed=1","--bic_var=dev",
      paste0("--outdir=", over, "/factorization_results/"),
      "--only_run=GLEANER_glmnet", "--K=5",paste0("--C=", over, "/sim1.c_matrix.txt"),
      "--no_plots")
  
  
  
  over="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/1b_overlap//V102_U102_MAF-mix_eur_N-50000_RHO-1b_high_mixed_p_No-1b_high_no/"
  t=c(paste0("--se_data=",over,"sim10.std_error.txt"),
      paste0("--beta_data=",over,"sim10.effect_sizes.txt"),
      "--seed=1","--bic_var=dev",
      paste0("--z_scores=",over,"sim10.z.txt"),
      paste0("--n_samples=",over,"sim10.N.txt"),
      paste0("--outdir=", over, "/factorization_results/sim10DEBUG"),
      "--only_run=GLEANER_glmnet", "--K=5",paste0("--C=", over, "/sim10.c_matrix.txt"))
  
  
  
#args <- parse_args(OptionParser(option_list=option_list),args=t)
args <- parse_args(OptionParser(option_list=option_list))

set.seed(args$seed)
#Read in the z-scores:
if(args$gwas_data != "")
{
	gwas_data <- args$gwas_data
	z_scores <- fread(gwas_data)
} else{
	beta <- fread(args$beta_data)
	se <- fread(args$se_data)
	#This could muff everything up....
	z_scores <- beta[,-1]/se[,-1]
	z_scores <- cbind(beta[,1], z_scores)
	colnames(z_scores) = c("ids", colnames(beta)[-1])
}

print(paste("Read in data of size", dim(z_scores)[1], "x",dim(z_scores)[2]))
old_rows <- dim(z_scores)[1]
z_scores <- z_scores %>% drop_na() #%>% separate(ids, into = c("chr", "pos"), sep = ":")
print(paste("Removed ", old_rows - dim(z_scores)[1], "SNPs containing missing data."))

#Which SNP got dropped?
#Import the trait names table: "trait map"
#Need to drop NAs too!
if(args$trait_names != "")
{
  trait_names <- scan(args$trait_names, what = character())
}else
{
  trait_names <- names(z_scores)[-1]
}

if(args$debug)
{
  z_scores <- z_scores[1:100,]
}
#Write out the id names for reference later on....

id_names <- z_scores %>% dplyr::select(ids)
#write_tsv(id_names, paste0(args$outdir, 'snp.ids.', args$seed, ".txt") )
write.table(id_names, paste0(args$outdir, 'snp.ids.', args$seed, ".txt"),  row.names = FALSE, quote = FALSE )
zstat_m <- as.matrix(z_scores %>% dplyr::select(-ids)) #sans the first columns, which are presumably snp id.
#zstat_m <- as.matrix(z_scores %>% select(-chr, -pos)) #sans the first columns, which are presumably snp id.

if(args$C != "")
{
  c.mat <- as.matrix(fread(args$C))
  message("Whitening")
}else{
  c.mat <- NULL
}

if(args$single_method != "")
{
  runs <- list()
  runs[[1]] <- runSingle(args$single_method, zstat_m, args$K, se_m = NULL, covar = c.mat)
  desc <- c(args$single_method)
  plotMatrices(runs[[1]], desc[1], trait_names, unlist(id_names))
  writeMatrices(runs[[1]], desc[1])
  writePerformance(runs[[1]], desc[1], zstat_m)
  quit()
}

#This is default. Run all the factorization methods we are considering.
if(args$all_run || args$only_run != "")
{

  beta <- fread(args$beta_data) %>% drop_na()
  se <- fread(args$se_data) %>% drop_na()

  if(args$debug)
  {
    beta <- beta[1:100,]
    se <- se[1:100,]
  }

    #TODO: fix this argument...
  #if(!args$no_column_id)
  if(any(grepl(pattern="rs",x=beta[,1])) | any(grepl(pattern="SNP",x=beta[,1])))
  {
    message("Removing 1st column, please verify this is correct")
    beta_m <- as.matrix(beta[,-1]) %>% apply(., 2, as.numeric)
    se_m <- as.matrix(se[,-1]) %>% apply(., 2, as.numeric)
    if(all(se[,1] != z_scores[,1]) || all(beta[,1] != z_scores[,1]))
    {
      print("Error in aligning")
      print("You likely don't have SNP identifiers on the rows")
      quit()
    }
  }else
  {
    beta_m <- as.matrix(beta)
    se_m <- as.matrix(se)
  }


  #check that the order is right
  runs <- list()
  #Run double PMA:
  K = args$K
  #Some methods take Z-scores directly, whearas others do not
  unscaled.methods <- c("PMA2", "backfit_noscale", "flashColumnVar","factorGo", "GLEANER","GLEANER_glmnet","GLEANER_glmnet_noCovar", "SVD_beta")
  scaled.methods <- c("sSVD", "PCA","SVD", "backfit","ssvd", "sPCA", "PCA_chooseK", "SVD_whiten")
  desc <- c("PMA2", "flashUbiqFilter", "flashUbiq_noscale", "ashr", "backfit_noscale", "backfit", "ssvd", "sPCA", "PCA", "flashColumnVar", "SVD_whiten")

  if(args$only_run != "")
  {
    desc <- str_split(args$only_run, ",")[[1]]
    message("running methods only specified,")
    print(paste(desc))
  }
  time.perf <- list()
  for(i in 1:length(desc))
  {
    meth = desc[i]
    message(paste0("Running ", meth))
    se = se_m
    effect.matrix <- beta_m
    if(meth %in% scaled.methods) #that is, it is a scaled method
    {
      #NOTE: this was a major change from before, I think the logic was wrong here.....
      se = NULL
      effect.matrix <- zstat_m
    }
    start_time <- Sys.time()
    if(grepl("gwasMF", meth) | grepl("GLEANER", meth))
    {
      #OCT 28- temporary hack to just run flash, trying to make a quick tweak....
      message("Skip")
      next;
      start_time <- Sys.time()
      message("Changing to converge at just 0.005")
      message(args$WLgamma)
	    run.data <- runSingle(meth, effect.matrix, K, se_m=se, covar = c.mat, bic.var = args$bic_var,
                            init.mat = args$init_mat, is.sim = TRUE, save.path = paste0(args$outdir, meth), 
	                          scale.mats = args$step_scaling,shrinkWL=args$WLgamma,conv_objective=0.005,min_bic_search_iter=5) #change on july 11,for testing#not accounting for BIC typoe here...
	    end_time <- Sys.time()
    } else
    {
      if(grepl("FLASH_SE$", meth))
      {
        run.data <- runSingle(meth, effect.matrix, K, z_path=args$z_scores, n_path=args$n_samples,se_m=se, covar = c.mat,savepath = paste0(args$outdir, meth)) #llast 2 args for factorgo only
      }else
      {
        #OCT 28- temporary hack to just run flash, trying to make a quick tweak....
        message("Skip")
        next;
      }
    }
    end_time <- Sys.time()
    time.perf[[meth]] <- end_time - start_time
    message("Runtime, ", time.perf[[meth]])
    if(!args$no_plots)
    {

    plotMatrices(run.data, meth,trait_names, unlist(beta[,1]))
    }
    writeAllOutputs(run.data, meth, run.data$X)
  }
  message("writing out time stamps to..", paste0(args$outdir, "time.stamps.Rdata"))
  save(time.perf , file = paste0(args$outdir, "time.stamps.Rdata"))
}
  

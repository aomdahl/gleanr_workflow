################################################################################################################################
## Ashton Omdahl, November 2024

##
################################################################################################################################
pacman::p_load(magrittr, dplyr, ggplot2, data.table, ggsignif, devtools, optparse, gleanr)

option_list <- list(
make_option("--run_data", type = "character", help = "Path to the RData file containing the true V matrix"),
make_option(c("--gwas_effects"), type = 'character', help = "Specify the Z or B file, depending on specified weighting scheme. First column is ids of each variant, column names specify the trait"),
make_option(c("--uncertainty"), type = 'character', help = "Specify the path to the SE or other uncertainty file, depending on the weightin scheme.irst column is ids of each variant, column names specify the trait"),
make_option(c("--trait_names"), type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the order in the input tables.", default = ""),
make_option(c("--weighting_scheme"), type = 'character', help = "Specify either Z, B, B_SE, B_MAF", default = "B_SE"),
make_option(c("--covar_matrix"), type = 'character', help = "Path to LDSC estimates of covariance effect. No adjustment made if none provided", default = ""),
make_option(c("-c", "--converged_obj_change"), type = 'numeric', help = "Specify the objective percent change required to achieve convergence", default = 0.001),
make_option(c("--outdir"), type = "character", help = "Where to write stuff out to,"),
make_option(c("--fixed_first"), type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE),
make_option(c("-K", "--nfactors"), type = "character", help = "specify the number of factors. Options are a number, or MAX, KAISER, K-2, CG", default = "GRID"),
make_option(c("--bic_var"), type = 'character', help = "Specify the bic method to use. Options are [sklearn_eBIC,sklearn,Zou_eBIC,Zou,dev_eBIC,dev, std,NONE]", default = "sklearn_eBIC"),
make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud"),#These related to covariance matrix tweaks
make_option(c("-s", "--sample_sd"), type="character", default= "", help="File containing the standard deviation of SNPs; if provided, used to scale LDSC gcov terms."),
make_option(c("-b", "--block_covar"), type = "numeric", default= 0.2, help="Specify the degree to which a block structure is enforced, by cluster distance. Default is 0.2"),
make_option(c("--covar_se_matrix"), type = "character", default="", help="Path to covar se matrix, needed if using Strimmer gamma specification."),
make_option(c("--seed"), type = "integer", default=1, help="Specify a seed to initialize simulations from."),
make_option(c("--model_select"), type = "logical", default=FALSE,action ="store_true", help="Specify this if you only wish to do the model selection step."),
make_option(c("-g", "--WLgamma"), type="character", default= "0",
            help="Specify the extent of the WL shrinkage on the input covariance matrix.\nCan pass as a number (1 for no covariance adjustment, 0 for full adjustment), or to specify a method (either MLE or Strimmer)"),
make_option("--n_permutations", type = "integer", help = "Number of permutations to perform", default = 100),
make_option("--shuff_type", type = "character", help = "Elements to shuffle, either `rows` or `rows_by_col`.", default = "rows_by_col"),
make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per study. TODO: Also has the flexibility to expand"),
#Extra arguments not used here, just included for continuity with gleanr preprocessing:
make_option(c("--drop_phenos"), type = "character", help = "Specify phenotypes to exclude (useful if encountering issues with covariance adjustment)", default = "")
)



t=c("--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.se.tsv",
    "--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.beta.tsv",
"--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete.trait_list.tsv",
    "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//ldsr_results/panUKBB_complete/summary_data/gcov_int.tab.csv",
"--outdir=/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/permute_testing/",
"--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete/sample_sd_report.tsv","--fixed_first",
"--WLgamma=Strimmer","--covar_se_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv",
"--nfactors=134", "--bic_var=sklearn_eBIC", "--converged_obj_change=0.001", "--WLgamma=Strimmer", "--model_select",
"--run_data=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")


args_input <- parse_args(OptionParser(option_list=option_list))#, args = t)
#args_input <- parse_args(OptionParser(option_list=option_list), args = t)
args <-fillDefaultSettings(args_input)
writeRunReport(args)

# Set output directory
outdir <- normalizePath(args$outdir, mustWork = FALSE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)




#Establish the internal settings
option <- readInSettings(args)
#Load in the data
input.dat <- readInData(args)
X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names; W_c <- input.dat$W_c; C <- input.dat$C_block
option$C <- C
option$WLgamma <- args$WLgamma
reg.vect <- prepRegressionElements(X,W,W_c,option) #returns sparse matrices since that saves a ton of memory by doing this once up front


# Load the true V matrix
load(args$run_data)  # Assumes it loads an object named `ret`
V_true <- ret$V

#Transition over to the full run:
option$K <- ncol(ret$V)
option$alpha1 <- ret$autofit_alpha[1]#This is applied to U
option$lambda1 <- ret$autofit_lambda[1] #This is applied to V
option$iter <- 0.5 #only fit U

# Display message for verbosity level
userMessage(args$verbosity, "Outcome data weighted and transformed for regression", thresh = 0)

# Initialize permutation storage
shuff_results_U <- list()
shuff_results_V <- list()

# Set random seed for reproducibility
set.seed(args$seed)

# Main permutation loop
for (i in 1:args$n_permutations) {
  message(sprintf("Running permutation %d/%d", i, args$n_permutations))
  
  # Shuffle V matrix, 
  if(args$shuff_type == "rows_by_col") #here shuffling each column separate
  {
    shuffled_V <- apply(V_true, 2, function(x) sample(x, length(x), replace = FALSE))
  }else
  {
    shuffled_V <- V_true[sample(1:nrow(V_true), nrow(V_true), replace = FALSE),]
  }
  
  
  # Run optimization
  #Note- with argument "option$iter = 0.5", the fitting shorts out and just returns U.
  ret_permute <- gwasML_ALS_Routine(option, X, W, W_c, shuffled_V, reg.elements = reg.vect, no.prune = TRUE)
  
  # Store results
  shuff_results_U[[i]] <- ret_permute
  shuff_results_V[[i]] <- shuffled_V
  # Save intermediate results every 100 iterations
  if (i %% 100 == 0 || i == args$n_permutations) {
    save_file <- file.path(outdir, sprintf("permute_%d.RData", i))
    save(shuff_results_U,shuff_results_V, file = save_file)
    message(sprintf("Saved results at iteration %d to %s", i, save_file))
    
    # Clear stored results to save memory
    shuff_results_U <- list()
    shuff_results_V <- list()
  }
}

message("Permutation process completed.")


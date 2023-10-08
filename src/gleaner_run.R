library(devtools)
library(optparse)

option_list <- list(
make_option(c("--gwas_effects"), type = 'character', help = "Specify the Z or B file, depending on specified weighting scheme. First column is ids of each variant, column names specify the trait"),
make_option(c("--uncertainty"), type = 'character', help = "Specify the path to the SE or other uncertainty file, depending on the weightin scheme.irst column is ids of each variant, column names specify the trait"),
make_option(c("--lambda_gc"), type = 'character', help = "Specify the path to the genomic correction coefficients. If none provided, none used", default = ""),
make_option(c("--trait_names"), type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the order in the input tables.", default = ""),
make_option(c("--weighting_scheme"), type = 'character', help = "Specify either Z, B, B_SE, B_MAF", default = "B_SE"),
make_option(c("--covar_matrix"), type = 'character', help = "Path to LDSC estimates of covariance effect. No adjustment made if none provided", default = ""),
make_option(c("-c", "--converged_obj_change"), type = 'numeric', help = "Specify the objective percent change required to achieve convergence", default = 0.001),
make_option(c("-i", "--niter"), type = 'numeric', help = "Cap the number of iterations on the matrix learning step", default = 300),
make_option(c("--outdir"), type = "character", help = "Source file location"),
make_option(c("--drop_phenos"), type = "character", help = "Specify phenotypes to exclude (useful if encountering issues with covariance adjustment)", default = ""),
make_option(c("--fixed_first"), type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE),
make_option(c("--debug"), type = "logical", help = "if want debug run", action = "store_true", default = FALSE),
make_option(c("--overview_plots"), type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE),
make_option(c("-K", "--nfactors"), type = "character", help = "specify the number of factors. Options are a number, or MAX, KAISER, K-2, CG", default = "MAX"),
make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand"),
make_option(c("--bic_var"), type = 'character', help = "Specify the bic method to use. Options are [sklearn_eBIC,sklearn,Zou_eBIC,Zou,dev_eBIC,dev, std]", default = "sklearn_eBIC"),
make_option(c("-p", "--param_conv_criteria"), type = 'character', help = "Specify the convergene criteria for parameter selection", default = "BIC.change"),
make_option(c("-r", "--rg_ref"), type = 'character', help = "Specify a matrix of estimated genetic correlations to initialize V from", default = ""),
make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud"),
#These related to covariance matrix tweaks
make_option(c("-s", "--sample_sd"), type="character", default= "", help="File containing the standard deviation of SNPs; if provided, used to scale LDSC gcov terms."),
make_option(c("-b", "--block_covar"), type = "numeric", default= 0.2, help="Specify the degree to which a block structure is enforced, by cluster distance. Default is 0.2"),
make_option(c("-g", "--WLgamma"), type="numeric", default= 0, help="Specify the extent of the WL shrinkage- higher is better!")
)

#Tools for debugging
#7/19/2023
# std BIC errors
t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//ukbb_benchmark.beta.tsv",
  "--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//ukbb_benchmark.se.tsv",
  "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//pheno.names.txt",
  "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark/comparative_analysis//BIC-std_K-MAX.DEBUG", 
  "--fixed_first",  "--nfactors=MAX", "--bic_var=std","--verbosity=1")
# no data to access error
t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//ukbb_benchmark.beta.tsv",
     "--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//ukbb_benchmark.se.tsv",
     "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//pheno.names.txt",
     "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark/comparative_analysis//BIC-std_K-MAX.DEBUG", 
     "--fixed_first",  "--nfactors=MAX", "--bic_var=dev","--verbosity=1")


t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//ukbb_benchmark.beta.tsv",
     "--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//ukbb_benchmark.se.tsv",
     "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//pheno.names.txt",
     "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark/comparative_analysis//sklearn_max_covar_manual", 
     "--fixed_first",  "--nfactors=MAX", "--bic_var=sklearn","--verbosity=1", 
     "--covar=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark/summary_data/gcov_int.tab.csv",
     "--genomic_correction=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark/summary_data/genomic_correction_intercept.tsv",
     "--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark/sample_sd_report.tsv")

args_input <- parse_args(OptionParser(option_list=option_list))#, args = t)
#TODO: cut oout this middle-man, make a more efficient argument input vector.

setwd("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/"); load_all()
args <-fillDefaultSettings(args_input)

writeRunReport(args)

option <- readInSettings(args)
option$alpha1 <- NA
option$lambda1 <- NA
outdir <-args$outdir

#Read in the hyperparameters to explore
hp <- readInParamterSpace(args)
input.dat <- readInData(args)
X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names; W_c <- input.dat$W_c; C <- input.dat$C_block
#Run the bic thing...

bic.dat <- getBICMatricesGLMNET(opath,option,X,W,W_c, all_ids, names)
save(bic.dat,option, file = paste0(outdir, "_bic_dat.RData"))
if(is.na(bic.dat$K))
{
  message("No signal detected under current settings. Program will end")
  quit()
}
#Transition over to the full run:
option <- bic.dat$options
option$K <- bic.dat$K
option$alpha1 <- bic.dat$alpha
option$lambda1 <- bic.dat$lambda

ret <- gwasML_ALS_Routine(option, X, W, W_c, bic.dat$optimal.v, maxK=bic.dat$K) #I like this better
ret[["snp.ids"]] <- all_ids
ret[["trait.names"]] <- names
save(ret,option, file = paste0(outdir, "_final_dat.RData"))

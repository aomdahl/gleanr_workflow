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
make_option(c("--bic_var"), type = 'character', help = "Specify the bic method to use. Options are [sklearn_eBIC,sklearn,Zou_eBIC,Zou,dev_eBIC,dev, std,NONE]", default = "sklearn_eBIC"),
make_option(c("-p", "--param_conv_criteria"), type = 'character', help = "Specify the convergene criteria for parameter selection", default = "BIC.change"),
make_option(c("-r", "--rg_ref"), type = 'character', help = "Specify a matrix of estimated genetic correlations to initialize V from", default = ""),
make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud"),
#These related to covariance matrix tweaks
make_option(c("-s", "--sample_sd"), type="character", default= "", help="File containing the standard deviation of SNPs; if provided, used to scale LDSC gcov terms."),
make_option(c("-b", "--block_covar"), type = "numeric", default= 0.2, help="Specify the degree to which a block structure is enforced, by cluster distance. Default is 0.2"),
make_option(c("--covar_se_matrix"), type = "character", default="", help="Path to covar se matrix, needed if using Strimmer gamma specification."),
make_option(c("-g", "--WLgamma"), type="character", default= "0", 
            help="Specify the extent of the WL shrinkage on the input covariance matrix.\nCan pass as a number (1 for no covariance adjustment, 0 for full adjustment), or to specify a method (either MLE or Strimmer)")
)

#Tools for debugging
#7/19/2023
# std BIC error

t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//ukbb_benchmark.beta.tsv",
     "--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//ukbb_benchmark.se.tsv",
     "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark//pheno.names.txt",
     "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark/comparative_analysis//sklearn_max_covar_manual", 
     "--fixed_first",  "--nfactors=MAX", "--bic_var=sklearn","--verbosity=1", 
     "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark/summary_data/gcov_int.tab.csv",
     "--genomic_correction=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark/summary_data/genomic_correction_intercept.tsv",
     "--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark/sample_sd_report.tsv")

t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.beta.tsv",
     "--uncertainty=/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.se.tsv",
     "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/pheno.names.txt",
     "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_GWAS_h2-0.1_rg-0.9/GLEANER_WL_0.7/", 
     "--fixed_first",  "--nfactors=MAX", "--bic_var=sklearn","--verbosity=1", 
     "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/ukbb_GWAS_h2-0.1_rg-0.9/summary_data/gcov_int.tab.txt",
     "--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/sample_sd_report.tsv",
     "--WLgamma=0.7")

#
t=c("--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete.se.tsv",
    "--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete.beta.tsv",
"--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete.trait_list.tsv", 
    "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//ldsr_results/panUKBB_complete/summary_data/gcov_int.tab.SELECTED.csv",
"--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/panUKBB_complete/",
"--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//sample_sd_report.SELECTED.tsv","--fixed_first", 
"--WLgamma=Strimmer","--covar_se_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.SELECTED.csv",
"--nfactors=KAISER")
 

t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2/ukbb_benchmark_2_conservative.beta.tsv",
     "--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//ukbb_benchmark_2_conservative.se.tsv",
     "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//pheno.names.txt",
     "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/comparative_analysis//sklearn_max_covar_manual", 
     "--fixed_first",  "--nfactors=KAISER", "--bic_var=sklearn","--verbosity=1",
     "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/ukbb_benchmark_2/summary_data/gcov_int.tab.csv",
     "--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2/sample_sd_report.tsv")

#finngen version of above:
t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/finngen_benchmark_2_conservative.beta.tsv",
     "--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//finngen_benchmark_2_conservative.se.tsv",
     "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//pheno.names.txt",
     "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/comparative_analysis//sklearn_max_covar_manual", 
     "--fixed_first",  "--nfactors=KAISER", "--bic_var=sklearn","--verbosity=1",
     "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark_2/summary_data/gcov_int.tab.csv",
     "--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/sample_sd_report.tsv")
args_input <- parse_args(OptionParser(option_list=option_list))#, args = t)
#args_input <- parse_args(OptionParser(option_list=option_list), args = t)
#TODO: cut oout this middle-man, make a more efficient argument input vector.

setwd("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/"); load_all()
args <-fillDefaultSettings(args_input)

writeRunReport(args)

#Establish the internal settings
option <- readInSettings(args)
option$alpha1 <- NA
option$lambda1 <- NA
outdir <-args$outdir

#Read in the hyperparameters to explore
hp <- readInParamterSpace(args)
#args$WLgamma <- 0
#args$block_covar <- 0.2 This is the issue
input.dat <- readInData(args)
X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names; W_c <- input.dat$W_c; C <- input.dat$C_block
option$C <- C
option$WLgamma <- args$WLgamma

#Prep regression datasets (weights, x, etc.)
#Keep these things in memory so not recreating each time
reg.vect <- prepRegressionElements(X,W,W_c) #returns sparse matrices since that saves a ton of memory.
#Loading that all in each time was slow...
print("Data prepped for regression")
lobstr::mem_used()
#What I should do is initialize a GLEANER object, and have all these things in there....
#Can I make the regvect elements sparse?
#Perform sparsity selection,
bic.dat <- getBICMatricesGLMNET(opath,option,X,W,W_c, all_ids, names,reg.elements=reg.vect)
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

#Perform optimization
ret <- gwasML_ALS_Routine(option, X, W, W_c, bic.dat$optimal.v, maxK=bic.dat$K,reg.elements=reg.vect) #I like this better
ret[["snp.ids"]] <- all_ids
ret[["trait.names"]] <- names
save(ret,option, file = paste0(outdir, "_final_dat.RData"))

#Write output files
write.table(ret$V, file = paste0(outdir, "latent.factors.txt"), quote = FALSE, row.names = ret$trait.names, sep = " ")
write.table(ret$U, file = paste0(outdir, "latent.loadings.txt"), quote = FALSE, row.names = ret$snp.ids, sep = " ")

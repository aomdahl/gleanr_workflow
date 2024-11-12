################################################################################################################################
## Ashton Omdahl, November 2024
## Wrapper script to run gleanr fro mthe command line
## Input:
##              - X: the matrix to be decomposed (N x M, N is the number of SNPs, M is the number of traits/studies)
##              - V: learned factor matrix  (M x K, M is the number of traits, K is the number of factors)
##              - W: the weight matrix, same size as in X
##              - option: a list with parameters including lambda1, the l1 penalty parameter for the factor matrix (F)
##
## Output:
##              - A text file containing matrix factors (`latent.factors.txt`) and latent loadings (`latent.loadings.txt`)
##              - Three R data files corresponding to sparsity parameter selection (`*_bic_dat.RData`), K_init selection (`*XXX`), and model fitting (`*_final_dat.RData`)
##
## Example of usage:
## Rscript src/gleaner_run.R --gwas_effects target_studies.beta.tsv --uncertainty target_studies.se.tsv --trait_names pheno.names.txt --outdir output_dir/prefix \
## --fixed_first -K "GRID" -v 1 \
## --covar_matrix target_studies_gcov_int.tab.csv \
## --sample_sd target_studies_sample_sd_report.tsv  \
## --WLgamma Strimmer \
## --covar_se_matrix target_studies_gcov_int_se.tab.csv
##
################################################################################################################################



library(devtools)
library(optparse)

writeOutResults <- function(V,U,traits, snps,outdir)
{
  if(length(ncol(V)) == 0) #If we get nothing back
  {
    file.create( paste0(outdir, "latent.factors.txt"))
    file.create(paste0(outdir, "latent.loadings.txt"))
  }else
  {
    V.df <- data.frame(traits, V) %>% magrittr::set_colnames(c("Study", paste0("V", 1:ncol(V))))
    write.table(V.df, file = paste0(outdir, "latent.factors.txt"), quote = FALSE, row.names = FALSE, sep = " ")
    
    U.df <- data.frame(snps, U) %>% magrittr::set_colnames(c("SNP", paste0("U", 1:ncol(V))))
    write.table(U.df, file = paste0(outdir, "latent.loadings.txt"), quote = FALSE, row.names = FALSE, sep = " ")
    
    write.table(traits, file =paste0(outdir, "trait_out_order.txt"), quote = FALSE, row.names = FALSE,col.names = FALSE)
  }

}




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
make_option(c("-K", "--nfactors"), type = "character", help = "specify the number of factors. Options are a number, or MAX, KAISER, K-2, CG", default = "GRID"),
make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand"),
make_option(c("--bic_var"), type = 'character', help = "Specify the bic method to use. Options are [sklearn_eBIC,sklearn,Zou_eBIC,Zou,dev_eBIC,dev, std,NONE]", default = "sklearn_eBIC"),
make_option(c("-p", "--param_conv_criteria"), type = 'character', help = "Specify the convergene criteria for parameter selection", default = "BIC.change"),
make_option(c("-r", "--rg_ref"), type = 'character', help = "Specify a matrix of estimated genetic correlations to initialize V from", default = ""),
make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud"),
make_option(c("-n", "--ncores"), type="integer", default= 1, help="Do you want to run on multiple cores? Only influences the K initialization step at this point."),
#These related to covariance matrix tweaks
make_option(c("-s", "--sample_sd"), type="character", default= "", help="File containing the standard deviation of SNPs; if provided, used to scale LDSC gcov terms."),
make_option(c("-b", "--block_covar"), type = "numeric", default= 0.2, help="Specify the degree to which a block structure is enforced, by cluster distance. Default is 0.2"),
make_option(c("--covar_se_matrix"), type = "character", default="", help="Path to covar se matrix, needed if using Strimmer gamma specification."),
make_option(c("--intermediate_bic"), type = "character", default="", help="OPTIONAL:Path to intermediate BIC file, will initiate factorization from here. "),
make_option(c("--subset_seed"), type = "numeric", default=-1, help="Specify a seed to subset the data for a variant test. For debugging purposes."),
make_option(c("--model_select"), type = "logical", default=FALSE,action ="store_true", help="Specify this if you only wish to do the model selection step."),
make_option(c("-g", "--WLgamma"), type="character", default= "0", 
            help="Specify the extent of the WL shrinkage on the input covariance matrix.\nCan pass as a number (1 for no covariance adjustment, 0 for full adjustment), or to specify a method (either MLE or Strimmer)")
)


#Tools for debugging
#7/19/2023
# std BIC error
#
t=c("--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.se.tsv",
    "--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.beta.tsv",
"--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete.trait_list.tsv", 
    "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//ldsr_results/panUKBB_complete/summary_data/gcov_int.tab.csv",
"--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/panUKBB_complete_41K_sklearn_eBIC/DEBUG",
"--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete/sample_sd_report.tsv","--fixed_first", 
"--WLgamma=Strimmer","--covar_se_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv",
"--nfactors=34", "--bic_var=sklearn_eBIC", "--converged_obj_change=0.001", "--WLgamma=Strimmer", "--model_select")

t=c("--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete_clumped_r2-0.2.se.tsv",
    "--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete_clumped_r2-0.2.beta.tsv",
    "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//panUKBB_complete.trait_list.tsv", 
    "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//ldsr_results/panUKBB_complete/summary_data/gcov_int.tab.csv",
    "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/panUKBB_complete_41K_sklearn_eBIC/DEBUG",
    "--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//gwas_extracts/panUKBB_complete//sample_sd_report.tsv","--fixed_first", 
    "--WLgamma=Strimmer","--covar_se_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv",
    "--nfactors=134","--bic_var=sklearn_eBIC")

t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2/ukbb_benchmark_2_conservative.beta.tsv",
     "--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//ukbb_benchmark_2_conservative.se.tsv",
     "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2//pheno.names.txt",
     "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/DEBUG//sklearn_max_covar_manual", 
     "--fixed_first",  "--nfactors=GRID", "--bic_var=sklearn","--verbosity=1",
     "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/ukbb_benchmark_2/summary_data/gcov_int.tab.csv",
     "--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2/sample_sd_report.tsv")

#finngen version of above:
t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/finngen_benchmark_2_conservative.beta.tsv",
     "--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//finngen_benchmark_2_conservative.se.tsv",
     "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//pheno.names.txt",
     "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/dev_max_covar_manual", 
     "--fixed_first",  "--nfactors=GRID", "--bic_var=sklearn_eBIC","--verbosity=1","--ncores=1",
     "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark_2/summary_data/gcov_int.tab.csv",
     "--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/sample_sd_report.tsv",
     "--WLgamma=Strimmer", "--covar_se_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark_2/summary_data/gcov_int_se.tab.csv")
args_input <- parse_args(OptionParser(option_list=option_list))#, args = t)
#args_input <- parse_args(OptionParser(option_list=option_list), args = t)
#TODO: cut oout this middle-man, make a more efficient argument input vector.

#setwd("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/"); load_all()
library(gleanr)
args <-fillDefaultSettings(args_input)

writeRunReport(args)

#Establish the internal settings
option <- readInSettings(args)
option$alpha1 <- NA
option$lambda1 <- NA
outdir <-args$outdir
#option$std_y <- TRUE
#Read in the hyperparameters to explore
hp <- readInParamterSpace(args)
#args$WLgamma <- 0
#args$block_covar <- 0.2 This is the issue
input.dat <- readInData(args)
X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names; W_c <- input.dat$W_c; C <- input.dat$C_block
option$C <- C
option$WLgamma <- args$WLgamma

if(args$subset_seed > -1)
{
  message("Subsetting randomly, seed ", args$subset_seed)
  set.seed(args$subset_seed)
  sub.samples <- sample(1:nrow(X),args$subset_seed) #the seed is both thte seed and the number of samples
  X <- X[sub.samples,]
  W <- W[sub.samples,]
  all_ids <- all_ids[sub.samples]
}


#Prep regression datasets (weights, x, etc.)
#Keep these things in memory so not recreating each time
reg.vect <- prepRegressionElements(X,W,W_c,option) #returns sparse matrices since that saves a ton of memory by doing this once up front
#Loading that all in each time was slow...
print("Data prepped for regression")
lobstr::mem_used()
opath=outdir #need to get rid of this stupidity, verify its okay.
#What I should do is initialize a GLEANER object, and have all these things in there....
#Can I make the regvect elements sparse?
#Perform sparsity selection,

if(args$model_select)
{
	message("only doing model selection")
	print(option$K)
}

if(args$debug)
{
  save(args, option, outdir, X,W,all_ids,names,W_c,C,reg.vect,option, file = paste0(outdir, "_debug_data.RData"))
  quit()
}

if(args$intermediate_bic == "")
{
  bic.dat <- getBICMatricesGLMNET(outdir,option,X,W,W_c, all_ids, names,reg.elements=reg.vect)
  save(bic.dat,option, file = paste0(outdir, "_bic_dat.RData"))
}else
{
  message("Using initialized version..")
  #Need to update some of the options, if there are any differences....
  load(args_input$intermediate_bic)
  option$conv0 <- args$converged_obj_change
  message("Starting at ", bic.dat$K)
}

if(args$model_select)
{
	message("Model selection complete, terminating")
	quit()
}
if(is.na(bic.dat$K))
{
  message("No signal detected under current settings. Program will end")
  writeOutResults(bic.dat$optimal.v, bic.dat$optimal.u, names, all_ids, outdir)
  ret <- list("V"=bic.dat$optimal.v, "U"=bic.dat$optimal.u,"snp.ids"=all_ids, "trait.names"=names)
  save(ret,option, file = paste0(outdir, "_final_dat.RData"))
  
  
  quit()
}



#Transition over to the full run:
message("Convergence set to ",  option$conv0 )
option$K <- bic.dat$K
option$alpha1 <- bic.dat$alpha #This is applied to U
option$lambda1 <- bic.dat$lambda #This is applied to V

#Perform optimization
ret <- gwasML_ALS_Routine(option, X, W, W_c, bic.dat$optimal.v, maxK=bic.dat$K,reg.elements=reg.vect) #I like this better
#Next, try it with them swapped......
#Testing init with u: ret.u <- gwasML_ALS_Routine(option, X, W, W_c, bic.dat$optimal.v, maxK=bic.dat$K,reg.elements=reg.vect, preU=bic.dat$optimal.u)
ret[["snp.ids"]] <- all_ids
ret[["trait.names"]] <- names
save(ret,option, file = paste0(outdir, "_final_dat.RData"))


writeOutResults(ret$V,ret$U,ret$trait.names, ret$snp.ids,outdir)

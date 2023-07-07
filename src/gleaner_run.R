library(devtools)
setwd("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/"); load_all()


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
make_option(c("--drop_phenos"), type = "character", help = "Specify phenotypes to exclude (useful if encountering issues with covariance adjustment)"),
make_option(c("--fixed_first"), type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE),
make_option(c("--debug"), type = "logical", help = "if want debug run", action = "store_true", default = FALSE),
make_option(c("--overview_plots"), type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE),
make_option(c("-K", "--nfactors"), type = "character", help = "specify the number of factors", default = "MAX"),
make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand"),
make_option(c("--bic_var"), type = 'character', help = "Specify the bic method to use. Options are [sklearn_eBIC,sklearn,Zou_eBIC,Zou,dev_eBIC,dev, std]", default = "sklearn_eBIC"),
make_option(c("-pc", "--param_conv_criteria"), type = 'character', help = "Specify the convergene criteria for parameter selection", default = "BIC.change"),
make_option(c("-rg", "--rg_ref"), type = 'character', help = "Specify a matrix of estimated genetic correlations to initialize V from", default = ""),
make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud")
)
args_input <- parse_args(OptionParser(option_list=option_list))#, args = t)
#TODO: cut oout this middle-man, make a more efficient argument input vector.
args <-fillDefaultSettings(args_input)

writeRunReport(args)

option <- readInSettings(args)
option$alpha1 <- NA
option$lambda1 <- NA
outdir <-args$outdir

#Read in the hyperparameters to explore
hp <- readInParamterSpace(args)
input.dat <- readInData(args)
X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names; W_c <- input.dat$W_c
#Run the bic thing...

bic.dat <- getBICMatricesGLMNET(opath,option,X,W,W_c, all_ids, names)
save(bic.dat,option, file = paste0(outdir, "_bic_dat.RData"))

#Transition over to the full run:
option <- bic.dat$options
option$K <- bic.dat$K
option$alpha1 <- bic.dat$alpha
option$lambda1 <- bic.dat$lambda

ret <- gwasML_ALS_Routine(option, X, W, W_c, bic.dat$optimal.v, maxK=bic.dat$K) #I like this better
ret[["snp.ids"]] <- all_ids
ret[["trait.names"]] <- names
save(ret.dat,option, file = paste0(outdir, "/final_dat.RData"))

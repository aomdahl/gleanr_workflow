#Source everything you need:
#pacman::p_load(tidyr, plyr, dplyr, ggplot2, stringr, penalized, cowplot, parallel, doParallel, logr, coop, data.table, glmnet, svMisc, nFactors)
#renv::init("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/renv_f")
message("Warning: all files will be put into the renv directory.")
suppressPackageStartupMessages(library("tidyr"))
#suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("penalized"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("doParallel"))
suppressPackageStartupMessages(library("logr"))
suppressPackageStartupMessages(library("coop"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("glmnet"))
suppressPackageStartupMessages(library("svMisc"))
suppressPackageStartupMessages(library("nFactors"))
#suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("optparse"))
dir ="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/"
#dir = "/Users/aomdahl/Documents/Research/LocalData/snp_networks/gwas_spMF/src/"
source(paste0(dir, "fit_F.R"))
source(paste0(dir, "update_FL.R"))
source(paste0(dir, "fit_L.R"))
source(paste0(dir, "plot_functions.R"))
source(paste0(dir, 'compute_obj.R'))
source(paste0(dir, 'buildFactorMatrices.R'))
source(paste0(dir, 'sparsity_scaler.R'))
source(paste0(dir, 'cophenetic_calc.R'))
source(paste0(dir, 'read_in_tools.R'))
source(paste0(dir, 'regressionUtils.R'))

quickSort <- function(tab, col = 1)
{
  tab[order(tab[,..col], decreasing = TRUE),]
}

updateStatement  <- function(l,a,l_og, a_og, run,time)
{

  log_print(paste0("Run complete. (Time: ", time, " min)"))
  log_print(paste0("Sparsity params: F ", l, ", L ", a))
  log_print(paste0("Sparsity scale: F ", l_og, ", L ", a_og))
  log_print(paste0("Current F sparsity: ", run$V_sparsity), console = FALSE)
  log_print(paste0("Current L sparsity: ", run$U_sparsity), console = FALSE)
  log_print(paste0("Number of active factors: ", ncol(run$V)))
}
GLOBAL_BASE <- 0
#For debug purposes, was running into issues with this on parallelied setups
reportOpenFiles <- function(option)
  {
    if(!option$debug)
    {
      return()
    }
    t = system("lsof -u aomdahl1", TRUE)
    message(paste0("Number of system loaded files now is..", length(t)))
    if(GLOBAL_BASE == 0)
    {
      GLOBAL_BASE <<- length(t)
    }
    if(length(t) > GLOBAL_BASE)
    {
      message('new added written out to:')
      #Create file with all the lines....
      debug.path = paste0(option$out, gsub(Sys.time(), pattern = " ", replacement="_"), ".txt")
      x <- lapply(t, function(x) gsub(x,pattern = "\\s+", replacement = " ")[[1]])
      y <- lapply(x, function(x) str_split(x[[1]], pattern = " ")[[1]])
      final <- data.frame(do.call("rbind", y))
      colnames(final) <- final[1,]
      final <- final[-1,]
      write.table(final, file= debug.path , quote = FALSE, row.names = FALSE)
      #for(i in (GLOBAL_BASE+1):length(t))
      #{
      #  print(t[[i]])
      #}

      GLOBAL_BASE <<- length(t)
    }

  }
option_list <- list(
make_option(c("--gwas_effects"), type = 'character', help = "Specify the Z or B file, depending on specified weighting scheme. First column is ids of each variant, column names specify the trait"),
make_option(c("--uncertainty"), type = 'character', help = "Specify the path to the SE or other uncertainty file, depending on the weightin scheme.irst column is ids of each variant, column names specify the trait"),
make_option(c("--lambda_gc"), type = 'character', help = "Specify the path to the genomic correction coefficients. If none provided, none used", default = ""),
make_option(c("--trait_names"), type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the orderr in the input tables.", default = ""),
make_option(c("--weighting_scheme"), type = 'character', help = "Specify either Z, B, B_SE, B_MAF", default = "B_SE"),
make_option(c("--init_F"), type = 'character', default = "ones_eigenvect", help = "Specify how the F matrix should be initialized. Options are [ones_eigenvect (svd(cor(z))_1), ones_plain (alls 1s), plieotropy (svd(cor(|z|))_1)]"),
make_option(c("--init_L"), type = 'character', help = "Specify this option to start by initializing L rather than F. Options are [random, pleiotropy]. Use empty string for none (default)", default = ""),
make_option(c("--alphas"), type = 'character', default = "", help = "Specify which alphas to do, all in quotes, separated by ',' character"),
make_option(c("--lambdas"), type = 'character', default = "", help = "Specify which lambdas to do, all in quotes, separated by ',' character"),
make_option(c("--scaled_sparsity"), type = "logical", action = "store_true", default = FALSE, help = "Specify this to scale the sparsity params by dataset to be between 0 and 1"),
make_option(c("--output"), type = "character", help = "Source file location"),
make_option(c("--calibrate_k"), type = "logical", help = "Give just an estimate of K", action = "store_true", default = FALSE),
make_option(c("--cores"), type = "integer", help = "Number of cores", default = 1),
make_option(c("--fixed_first"), type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE),
make_option(c("--debug"), type = "logical", help = "if want debug run", action = "store_true", default = FALSE),
make_option(c("--overview_plots"), type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE),
make_option(c("-k", "--nfactors"), type = "character", help = "specify the number of factors", default = "0"),
make_option(c("-i", "--niter"), type = "integer", help = "specify the number of iterations", default = 30),
make_option(c("--posF"), type = "logical", default = FALSE,  help = "Specify if you want to use the smei-nonnegative setup.", action = "store_true"),
make_option(c("--scale_n"), type = "character", default = "",  help = "Specify the path to a matrix of sample sizes if you want to scale by sqrtN as well as W"),
make_option(c("--IRNT"), type = "logical", default = FALSE,  help = "If you want to transform the effect sizes ", action = "store_true"),
make_option(c("--MAP_autofit"), type = "integer", default = -1,  help = "Specify if you want to the MAP autofit sparsity parameter to learn lambda, alpha"),
make_option(c("-f", "--converged_F_change"), type="double", default=0.02,help="Change in factor matrix to call convergence"),
make_option(c("-o", "--converged_obj_change"), type="double", default=1e-4,help="Relative change in the objective function to call convergence"),
make_option(c("--no_SNP_col"), type="logical", default= FALSE, action = "store_true", help="Specify this option if there is NOT first column there..."),
make_option(c("--parameter_optimization"), type="integer", default= 1, action = "store_true", help="Specify how many iterations to run to get the cophenetic optimization, etc. "),
make_option(c("--regression_method"), type="character", default= "penalized", help="Specify which regression method to use. Options are [penalized (Default), glmnet]"),
make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand"),
make_option(c("--epsilon"), type="double", default= 1e-8, help="The convergence criteria for the L1. If exploring the space, try making this larger to speed up runtime. "),
make_option(c("--subsample"), type="integer", default= 0, help="Specify if you wish to subsample from the full dataset, and if so by how much. For repeated runs, this will choose a different subsample each time."),
make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud"),
make_option(c("--auto_grid_search"), type="logical", default=FALSE, action = "store_true", help="Specify if you want to do an auatomatic grid search"),
make_option(c("-c","--covar_matrix"), type = "character", default = "", help = "Path to a covariance matrix. Should be pre-filtered to include only significant elements")
)
args <- parse_args(OptionParser(option_list=option_list))
message("Please make sure the first column of input data is SNP/RSIDs.")
if(args$MAP_autofit >-1 & args$auto_grid_search)
{
  message("Incompatible arguments (grid search and MAP-based autofit). Please try again.")
}

#lf <- log_open(paste0(args$output, "gwasMF_log.", Sys.Date(), ".txt"))
#TODO:
    #add in functionality in the event the objective begins to increase again. Don't want that....
  #finesse functionality to allow for trait-specific variance to be learned (?)
#easy enough to add in.
if(FALSE) #For debug functionality on MARCC- this is currently loading the udler data.
{
  datdir <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/"
  datdir <- "/Users/aomdahl/Documents/Research/LocalData/snp_networks/datasets/"
  args <- list()
  #args$gwas_effects <- paste0(datdir, "beta_signed_matrix.tsv")
  #args$uncertainty <-  paste0(datdir, "se_matrix.tsv")

  args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/B.tsv"
  args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/SE.tsv"
  args$fixed_first <- TRUE
  args$genomic_correction <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/Lambda_gc.tsv"
  args$nfactors <- 5
  args$calibrate_k <- FALSE
  args$subsample <- -1
  args$trait_names = ""
  args$overview_plots <- FALSE
  args$parameter_optimization <- 5
  args$niter <- 15
  args$covar_matrix <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/ukbb_GWAS_h2-0.1_rg-0.9/summary_data/gcov_int.tab.new_names.txt"
  #args$alphas <- "0.01,0.05,0.1,0.15,0.2" #using 0.00001 since it doesn't like 0
  #args$lambdas <- "0.01,0.05,0.1,0.15,0.2"
  args$alphas <- "1000" #using 0.00001 since it doesn't like 0
  args$lambdas <- "36"
  args$MAP_autofit <- -1
  args$auto_grid_search <- FALSE
    args$cores <- 1
  args$IRNT <- FALSE
  args$fixed_first <- TRUE
  args$weighting_scheme = "B_SE"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/scratch/testing_gwasMF_code/gls/"
  args$converged_obj_change <- 0.1
  args$scaled_sparsity <- FALSE
  args$posF <- FALSE
  #args$autofit <- 0
  args$init_F <- "ones_eigenvect"
  args$init_L <- ""
  args$epsilon <- 1e-8
  args$verbosity <- 1
  args$regression_method <- "penalized"
  args$scale_n = ""
  args$debug <- FALSE
  args$cores = 1
  args$IRNT <- FALSE

}
log.path <- paste0(args$output, "gwasMF_log.", Sys.Date(), ".txt")
print(log.path)
lf <- log_open(log.path, show_notes = FALSE)
options("logr.compact" = TRUE)
options("logr.notes" = FALSE)

output <- args$output

#Read in the hyperparameters to explore
hp <- readInParamterSpace(args)
if(args$MAP_autofit < 0 & !args$auto_grid_search)
{
  log_print("Input sparsity settings")
  log_print("Alphas: ")
  log_print(hp$a)
  log_print("Lambdas: ")
  log_print(hp$l)
}

#set up settings
option <- readInSettings(args)
#reportOpenFiles(option)
option$log <- log.path
log_print("Succesfully loaded program options....")
#reportOpenFiles(option)
#read in the data
input.dat <- readInData(args)
X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names; W_c <- input.dat$W_c
Sigma <- input.dat$C;
#make it a block matrix for this analysis...
if(is.null(Sigma))
{
	option$covar <- diag(ncol(X))
}else
{
	option$covar <- blockifyCovarianceMatrix(Sigma)
}
#Estimate sparsity space (doesn't work too well.)
max_sparsity <- approximateSparsity(X, W, option) #Has errors when you don't specify sparsity
if(option$calibrate_k)
{
  message("Writing out a sparsity file")
  fileConn<-file(paste0(output,"init_k.txt"))
  writeLines(as.character(max_sparsity$newK), fileConn)
  close(fileConn)
  quit()
}

#TODO: fix this- this is all nonsense. Metho doesn't work, maximum suggested is invalid. Need to include the weighting, update to new code base.
log_print("Estimating sparsity maximums based on an SVD approximation.")
log_print(paste0("Max alpha: ", max_sparsity$alpha))
log_print(paste0("Max lambda: ", max_sparsity$lambda))
option$K <- max_sparsity$newK #sets K if not provided.
message("THE ABOVE estimates are INCORRECT. NEED TO UPDATE CODE. TODO")
#reportOpenFiles(option)
#Use the max to set up autofit
if(option$MAP_autofit > -1)
{
  autofit.seed <- 0.001
  alphas <- c(max_sparsity$alpha * autofit.seed)
  lambdas <- c(max_sparsity$lambda * autofit.seed)
  option$max_lambda = max_sparsity$lambda
  option$max_alpha = max_sparsity$alpha
}

message(paste0("Number of SNPs for analysis: ", nrow(X)))
message(paste0("Number of traits for analysis: ", ncol(X)))

if(option$calibrate_sparsity)
{
  if(option$MAP_autofit > -1)
  {
    log_print("Incompatible 'MAP_autofit' and 'calibrate_sparsity' settings specified. Program will terminate")
    quit()
  }
  option$max_alpha <- max_sparsity$alpha
  option$max_lambda <- max_sparsity$lambda
  if(all(hp$a <= 1) & all(hp$a >= 0) & all(hp$l <= 1) & all(hp$l >= 0) & option$calibrate_sparsity)
  {
    alphas <- hp$a * max_sparsity$alpha
    lambdas <- hp$l * max_sparsity$lambda

    option$calibrate_sparsity <- FALSE
    updateLog("Scaled amounts are:", option)
    updateLog("    Lambdas:", option)
    updateLog(paste(round(lambdas, digits = 2)), option)
    updateLog("    Alphas:", option)
    updateLog(round(alphas, digits = 3), option)
  }else
  {
    message('Inputed sparsity parameters must be between 0 and 1 to use the "calibrate_sparsity option".')
    message("Program will terminate.")
    quit()
  }
} else{
  alphas <- hp$a
  lambdas <- hp$l
}

if(option$auto_grid_search)
{
  start <- Sys.time()
  message("initializing for automatic grid search")
  init.params <- GetCoarseSparsityParams(X,W,option, burn.in.reps = 5, logs = TRUE, n.points = 4)
  alphas <- init.params$alphas; lambdas <- init.params$lambdas
  message("Alphas")
  print(alphas)
  message("Lambdas")
  print(lambdas)
  time <-  round(difftime(Sys.time(), start, units = "mins"), digits = 3)
  log_print(paste0("Sparsity parameter estimation time: ", time))
}

opti=1
#Figure out the parallelism
#Basically, if splitting across all the cores results in less than 1000 SNPs/core, I should just split this loop across other cores, i.e.
#if 20 cores specified, and only 10,000 SNPs
#I should split this loop across 2 cores, and then gi

#args$parameter_optimization <- 5
for(opti in 1:args$parameter_optimization)
{
  message(paste0("On optimization run number: ", opti))
  if(option$subsample > 0)
  {
    message("We are subsampling to ", option$subsample)
    samp <- sample(1:nrow(X), option$subsample)
    X <- input.dat$X[samp,]; W <- input.dat$W[samp,]; all_ids <- input.dat$ids[samp]
  }


  #reportOpenFiles(option)
  #Store stats here
  type_ <- args$weighting_scheme
  run_stats <- list()
  name_list <- c()

  a_plot <- c()
  l_plot <- c()
  iter_count <- 1
  #In the case of autofit, you wouldn't go through a bunch of starting points, would you? maybe. Actually not a bad idea, huh.

  for(i in 1:length(alphas)){
    for (j in 1:length(lambdas)){
      #reportOpenFiles(option)
        a <- alphas[i]
        l <- lambdas[j]
        option[['alpha1']] <- as.numeric(a)
        option[['lambda1']] <- as.numeric(l)
      start <- Sys.time()
      run <- Update_FL(as.matrix(X), as.matrix(W), option)
      time <-  round(difftime(Sys.time(), start, units = "mins"), digits = 3)
      run_stats[[iter_count]] <- c(run, time)
      iter_count <- iter_count + 1
      #reportOpenFiles(option)
      #things for plots
      if(length(run) == 0)
      {
        fname = paste0(output, "K", option$K, "_A", a, "_L", l, "_","K", type_, ".NO_OUTPUT.png") #add in k
      } else
      {
        #reportOpenFiles(option)
        updateStatement(l,a,hp$l[j],hp$a[i], run,time)
        fname = paste0(output, "K", option$K,"_A", round(a, digits = 4), "_L", round(l, digits = 4), "_", type_,".", opti, ".png")
        title_ = paste0("A", round(a, digits = 4), "_L", round(l, digits = 4), " Type = ", args$weighting_scheme )
        p <- plotFactors(run$V,trait_names = names, title = title_, scale.cols = TRUE)
        ggsave(filename = fname, plot = p, device = "png", width = 10, height = 8)
        #Lazy bookeeping nonsense for plots
        name_list <- c(name_list,  paste0("K", option$K,"_A", round(a, digits = 6), "_L", round(l, digits = 6)))
        a_plot <- c(a_plot, a)
        l_plot <- c(l_plot, l)
        writeFactorMatrix(names, all_ids, run,  paste0("K", option$K,"_A", round(a, digits = 6), "_L", round(l, digits = 6), "_", type_, ".", opti), output)
      }

    }
  }
  #Save all the data...
  #We actually need output information
  save(run_stats, file = paste0(output, "runDat.", opti, ".RData"))
  log_print(paste0("Run data written out to: ", output, "runDat.RData"))
  #writeFactorMatrices(round(alphas, digits = 3), round(lambdas, digits = 3), names,all_ids, run_stats,output)
  check_stats = max(sapply(run_stats, length))
  reportOpenFiles(option)
  if(check_stats == 1)
  {
    log_print("No runs completed. Terminating")
    quit()
  }
  if(args$overview_plots)
  {
    valid_run_stats <- run_stats[which(sapply(run_stats, length) > 0)] #only do ones where a full run was completed.
    if(length(valid_run_stats) == 0)
    {
      log_print("Unable to write out images... no completed runs.")
      break
    }
    message('this code is still so buggy')
    factor_sparsities <- data.frame("alpha" = a_plot, "lambda" = l_plot,
                                    "sparsity" = unlist(lapply(valid_run_stats, function(x) x[[4]])),
                                    "runtime" = unlist(lapply(valid_run_stats, function(x) x[[7]])))
    loading_sparsities <- data.frame("alpha" = a_plot, "lambda" = l_plot,
                                     "sparsity" = unlist(lapply(valid_run_stats, function(x) x[[3]])),
                                     "runtime" = unlist(lapply(valid_run_stats, function(x) x[[7]])))
    message(paste0(output, "factor_sparsity", opti, ".png"))

    ggplot(data = factor_sparsities, aes(x = alpha, y= lambda, fill = sparsity)) + geom_tile() +
      theme_minimal(15) + scale_fill_gradient(low="white", high="navy") + ggtitle("Factor Sparsity")
    ggsave(paste0(output, "factor_sparsity.", opti, ".png"), device = "png")

    ggplot(data = loading_sparsities, aes(x = alpha, y= lambda, fill = sparsity)) + geom_tile() +
      theme_minimal(15) + scale_fill_gradient2(low="white", high="navy", limits = c(0,max(loading_sparsities$sparsity))) + ggtitle("Loading Sparsity")
    ggsave(paste0(output, "loading_sparsity", opti, ".png"), device = "png")

    #Plot the change in objective function for each one too....
   max_len <- max(unlist(lapply(valid_run_stats, function(x) length(x[[6]]))))
   res <- lapply(valid_run_stats, function(x) c(x[[6]], rep(NA,max_len- length(x[[6]]))))
  objectives <- data.frame("obj" = res)
   colnames(objectives) <- name_list
   objectives <- objectives %>% mutate("iteration" =1:nrow(objectives)) %>% reshape2::melt(., id.vars = "iteration")
   ggplot(data = objectives, aes(x = iteration, y = value, color = variable)) + geom_point() + geom_line() +
      theme_minimal(15) + ggtitle("Objective function") + ylab("Objective score")
    ggsave(paste0(output, "objective.", opti, ".png"), device = "png")
  }
  rm(run_stats)
  log_print(paste0("Completed stability iteration number ", opti))
  log_print("----------------------------------")
  log_print("")
}
####Propose the optimal solution
##TODO: this is not right. Have updated code that does this better, use that.
param.combos <- unlist(lapply(alphas, function(x)
  lapply(lambdas, function(y) paste0("A", x, "_L", y))))

original.combos <- unlist(lapply(hp$a, function(x)
  lapply(hp$l, function(y) paste0("A", x, "_L", y))))

if(args$parameter_optimization > 1)
{
  #Go through all the attempted parameter options
  fmats <- list()
  #lmats <- list()
  #initialize storage
  len.tries = length(alphas) * length(lambdas)
  for(i in 1:len.tries)
  {
    fmats[[i]] <- list()
    #lmats[[i]] <- list()
  }

  for(opti in 1:args$parameter_optimization)
  {
    #drop ones with only one factor
    dpath = paste0(output, "runDat.", opti, ".RData")
    load(dpath)
    for(paramcombo in 1:len.tries)
    {
      fmats[[paramcombo]][[opti]] <- run_stats[[paramcombo]][[1]]
      #lmats[[paramcombo]][[opti]] <- run_stats[[paramcombo]][[2]]
    }
    rm(run_stats)
  }
  coph.results <- lapply(fmats, copheneticAssessment) #give a list of the 30
  corr.results <- lapply(fmats, correlationAssessment) #give a list of the 30

  #plot cophenetic results
  top.result <- which.max(coph.results)
  message("Based on finished results, we recommend raw ALPHA, LAMBDA: ", param.combos[top.result])
  message("This corresponds to scaled values of ALPHA, LAMBDA: ", original.combos[top.result])

  #report the average across these runs for that one.
  out <- data.frame("Coph" = unlist(coph.results), "Corr_frob" = unlist(corr.results), "R" =param.combos, "R_og" =  original.combos) %>% arrange(Coph, Corr_frob)
  write.table(out, file =  paste0(output, "performance_summary_file.txt"), sep = "\t", quote = FALSE, row.names=FALSE)
  print(out %>% filter(Coph > 0.9) %>% arrange(Corr_frob))

}

log_print("Option settings",console = FALSE)
log_print(option,console = FALSE)
log_close()

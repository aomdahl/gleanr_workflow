#This version is to run functions from R directly.
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("plyr")) 
suppressPackageStartupMessages(library("dplyr")) 
suppressPackageStartupMessages(library("stringr")) 
suppressPackageStartupMessages(library("penalized")) 
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
  log_print(paste0("Current F sparsity: ", run$F_sparsity), console = FALSE)  
  log_print(paste0("Current L sparsity: ", run$L_sparsity), console = FALSE)  
  log_print(paste0("Number of active factors: ", ncol(run$F)))
}

#ACtual functions to run the thing
runOptimizingRuns <- function(num_iter, X, W, option)
{
    factor.matrices <- list()
    for(opti in 1:args$parameter_optimization) 
    {
        message(paste0("On optimization run number: ", opti))
        if(option$subsample > 0)
          {
            message("We are subsampling to ", option$subsample)
            samp <- sample(1:nrow(X), option$subsample)
            X.p <- X[samp,]; W.p <- W[samp,]; 
          } else{
            X.p <- X; W.p <- W
          }
        run_stats <- list()
        name_list <- c()
        a_plot <- c()
        l_plot <- c()
        iter_count <- 1
        for(i in 1:length(alphas)){
            for (j in 1:length(lambdas)){
                a <- alphas[i]
                l <- lambdas[j]
                option[['alpha1']] <- as.numeric(a)
                option[['lambda1']] <- as.numeric(l)
                start <- Sys.time()
                run <- Update_FL(as.matrix(X.p), as.matrix(W.p), option)
                time <-  round(difftime(Sys.time(), start, units = "mins"), digits = 3)
                run_stats[[iter_count]] <- c(run, time)
                iter_count <- iter_count + 1
            
            }

        }
            check_stats = max(sapply(run_stats, length))
            if(check_stats == 1)
            {
                log_print("No runs completed this round. Terminating")    
            } else{
                factor.matrices[[opti]] <- run_stats
            }
    }
    return(factor.matrices)
}

chooseBestFactorization <- function(num_iter,solutions, alphas, lambdas)
{
    #Go through all the attempted parameter options
    fmats <- list()
    #initialize storage
    len.tries = length(alphas) * length(lambdas)
    for(i in 1:len.tries)
    {
        fmats[[i]] <- list()
        #lmats[[i]] <- list()
    }
  
  for(opti in 1:num_iter)
  {
    #drop ones with only one factor
    run_stats <- solutions[[opti]]
    for(paramcombo in 1:len.tries)
    {
      fmats[[paramcombo]][[opti]] <- run_stats[[paramcombo]][[1]]
    }
  }
    coph.results <- lapply(fmats, copheneticAssessment) 
    corr.results <- lapply(fmats, correlationAssessment)
  
  #plot cophenetic results
  top.result <- which.max(coph.results)
    param.combos <- unlist(lapply(alphas, function(x) 
    lapply(lambdas, function(y) paste0("A", x, "_L", y))))

    original.combos <- unlist(lapply(hp$a, function(x) 
    lapply(hp$l, function(y) paste0("A", x, "_L", y))))
  #report the average across these runs for that one.
  out <- data.frame("Coph" = unlist(coph.results), "Corr_frob" = unlist(corr.results), "R" =param.combos, "R_og" =  original.combos) %>% 
    filter(Coph > 0.9) %>% arrange(Corr_frob)
    if(nrow(out) == 0)
    {
        message("No hits meet criteria. Selecting one with max Cophenetic correlation")
       out <- data.frame("Coph" = unlist(coph.results), "Corr_frob" = unlist(corr.results), "R" =param.combos, "R_og" =  original.combos)) %>%
       arrange(-Coph)
    }
        best <- str_split(out$R[1],pattern = "_")[[1]]
        ba <- gsub(t[1], pattern = "A", replacement = ""); bl <- gsub(t[1], pattern = "L", replacement = "") 

    return(list('alpha' = as.numeric(ba), 'lambda' = as.numeric(bl)))
}


GLOBAL_BASE <- 0
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
make_option(c("--cores"), type = "integer", help = "Number of cores", default = 1),
make_option(c("--fixed_first"), type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE),
make_option(c("--debug"), type = "logical", help = "if want debug run", action = "store_true", default = FALSE),
make_option(c("--overview_plots"), type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE),
make_option(c("-k", "--nfactors"), type = "character", help = "specify the number of factors", default = "15"),
make_option(c("-i", "--niter"), type = "integer", help = "specify the number of iterations", default = 30),
make_option(c("--posF"), type = "logical", default = FALSE,  help = "Specify if you want to use the smei-nonnegative setup.", action = "store_true"),
make_option(c("--scale_n"), type = "character", default = "",  help = "Specify the path to a matrix of sample sizes if you want to scale by sqrtN as well as W"),
make_option(c("--IRNT"), type = "logical", default = FALSE,  help = "If you want to transform the effect sizes ", action = "store_true"),
make_option(c("--MAP_autofit"), type = "integer", default = -1,  
            help = "Specify if you want to autofit sparsity parameter for the whole thing using the MAP approach. 1 is std (both lambda, alpha at once), 0 is adjusted, -1 is none (default", default = -1),
make_option(c("-f", "--converged_F_change"), type="double", default=0.02,help="Change in factor matrix to call convergence"),
make_option(c("-o", "--converged_obj_change"), type="double", default=1,help="Relative change in the objective function to call convergence"),
make_option(c("--no_SNP_col"), type="logical", default= FALSE, action = "store_true", help="Specify this option if there is NOT first column there..."),
make_option(c("--parameter_optimization"), type="integer", default= 1, action = "store_true", help="Specify how many iterations to run to get the cophenetic optimization, etc. "),
make_option(c("--regression_method"), type="character", default= "penalized", help="Specify which regression method to use. Options are [penalized (Default), glmnet]"),
make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand"),
make_option(c("--epsilon"), type="double", default= 1e-8, help="The convergence criteria for the L1. If exploring the space, try making this larger to speed up runtime. "),
make_option(c("--subsample"), type="integer", default= 0, help="Specify if you wish to subsample from the full dataset, and if so by how much. For repeated runs, this will choose a different subsample each time."),
make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud")
)
#this assumes, X, W, etc. all read in already
gwasMF <- function(X, W, snp_ids, trait_names, init_F = "ones_eigenvect", lambda_gc=1,weighting_scheme="B_SE", 
        alphas = c(0.001, 0.005, 0.01,0.05,0.1), lambdas = c(0.001, 0.005, 0.01,0.05,0.1), scaled_sparsity = TRUE, 
        cores = 1, fixed_first = FALSE, K=0, niter=15, scale_n= "", obj_change = 1, param_opti=10, verbose=FALSE)
        {
            args <- list()
            #Input data
            args$gwas_effects <-X; args$uncertainty <- W;  args$scale_n <- scale_n; args$genomic_correction <- lambda_gc
            args$fixed_first <- fixed_first; args$nfactors <- K
            args$trait_names = trait_names
            args$output <- ""      
            #Optimization settings
            args$niter <- niter; args$alphas <- alphas; args$lambdas <- lambdas; args$MAP_autofit <- FALSE; args$cores <- cores; args$scaled_sparsity <- TRUE;
            #convergence settings
            args$converged_obj_change <- obj_change; args$epsilon <- 1e-8
            #initialization settings
            args$init_F <- init_F; args$init_L <- ""; args$IRNT <- FALSE; args$weighting_scheme = weighting_scheme
            args$posF <- FALSE
            args$verbosity <- ifelse(verbose, 1,0)


            hp <- readInParamterSpace(args)
            option <- readInSettings(args)
            message("Succesfully loaded program options....")
            all_ids <- snp_ids
            names <- trait_names
            #Need to check and modify X, W, etc. according to the inputs we get.
            #TODO- fix W, X, etc.
            message("Setting up sparsity parameter grid...")
            max_sparsity <- approximateSparsity(X, W, option)
            option$K <- max_sparsity$newK #sets K if not provided.
            if(scaled_sparsity)
            {
                if(option$calibrate_sparsity)
                    {
                    option$max_alpha <- max_sparsity$alpha
                    option$max_lambda <- max_sparsity$lambda
                    if(all(hp$a <= 1) & all(hp$a >= 0) & all(hp$l <= 1) & all(hp$l >= 0) & option$calibrate_sparsity)
                    {
                        alphas <- hp$a * max_sparsity$alpha
                        lambdas <- hp$l * max_sparsity$lambda
                        
                        option$calibrate_sparsity <- FALSE
                    }else
                    {
                        message('Inputed sparsity parameters must be between 0 and 1 to use the "calibrate_sparsity option".')
                        quit()
                    }
                    } else{
                    alphas <- hp$a
                    lambdas <- hp$l
                    }

            }

        }
        quick.run <- runOptimizingRuns(param_opti, X, W, option)
        optimal.solution <- chooseBestFactorization(param_opti, quick.run,alphas, lambdas)
        final.run <- runOptimizingRuns(1, )


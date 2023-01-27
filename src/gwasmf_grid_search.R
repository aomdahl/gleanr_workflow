#Source everything you need:

SourcePackages <- function()
{
  dir ="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/"
  source(paste0(dir, "fit_F.R"))
  source(paste0(dir, "update_FL.R"))
  source(paste0(dir, "fit_L.R"))
  source(paste0(dir, "plot_sfunctions.R"))
  source(paste0(dir, 'compute_obj.R'))
  source(paste0(dir, 'buildFactorMatrices.R'))
  source(paste0(dir, 'sparsity_scaler.R'))
  source(paste0(dir, 'cophenetic_calc.R'))
  source(paste0(dir, 'read_in_tools.R'))
  source(paste0(dir, 'regressionUtils.R'))
  source(paste0(dir, 'pve.R'))
}

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


initializeGwasMF <- function(X,W, snp.ids, trait.names, K=0)
{
  pacman::p_load(tidyr, plyr, dplyr, ggplot2, stringr, 
                 penalized, cowplot, parallel, doParallel, logr, 
                 coop, data.table, glmnet, nFactors)
  
  dir ="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/"
  source(paste0(dir, "gwasmf_grid_search.R"))
  SourcePackages()
  
  args <- defaultSettings(K=K)
  option <- readInSettings(args)
  option$swap <- FALSE
  option$alpha1 <- 1e-10
  option$lambda1 <- 1e-10
  output <- args$output
  log.path <- paste0(args$output, "gwasMF_log.", Sys.Date(), ".txt")
  lf <- log_open(log.path, show_notes = FALSE)
  options("logr.compact" = TRUE)
  options("logr.notes" = FALSE)
  #Read in the hyperparameters to explore
  hp <- readInParamterSpace(args)
  all_ids <-snp.ids; names <- trait.names
  initk <- option$K
  #1/16: testing out a new idea on the sims. Might be too much, but whatever.
  if(initk == 0)
  {
  message('Initializing X to the max -1')
  option$K <- ncol(X)-1
  }
  #Run the bic thing...
  option$V <- FALSE
  option$fixed_ubiqs <- TRUE
  list("args" = args, "options" = option, 
       "hp" = hp, "all_ids" = snp.ids, "names" = trait.names, "initk" = initk)
}

#getGridMatrices(args, option, alphas, lambdas, param.iters)
getGridMatrices <- function(args, option, alphas, lambdas, param.iters)
{
  args$parameter_optimization <- param.iters
  all.runs <- list()
  run.ids <- list()
  #Change on 1/27- put the replicates within the loop, so its easier on the downstream end...

    #In the case of autofit, you wouldn't go through a bunch of starting points, would you? maybe. Actually not a bad idea, huh.
    for(i in 1:length(alphas)){
      for (j in 1:length(lambdas)){
        a <- alphas[i]
        l <- lambdas[j]
        option[['alpha1']] <- as.numeric(a)
        option[['lambda1']] <- as.numeric(l)
        run_stats <- list()
        for(opti in 1:args$parameter_optimization) 
        {
          message(paste0("On optimization run number: ", opti))
          if(option$subsample > 0)
          {
            message("We are subsampling to ", option$subsample)
            samp <- sample(1:nrow(X), option$subsample)
            X <- input.dat$X[samp,]; W <- input.dat$W[samp,]; all_ids <- input.dat$ids[samp]
            #This works, but need to adjust the parameters. COuld implement okay, but not doing now.f
          }
        run_stats[[opti]] <- Update_FL(as.matrix(X), as.matrix(W), option)
        }
        run.id <- paste0(option$K, "_", a, "_", l)
        run.ids <- c(run.ids, run.id)
        save(run_stats, file = paste0(args$output, run.id, ".RData"))
        all.runs[[run.id]] <- run_stats
    }
    
  }
    return(all.runs)
}


#helper function to run the whole thing 
#gwasMFBIC(X,W, paste0("rs", 1:nrow(X)), paste0("T",1:ncol(X)), K=K,...)
#runSingle(meth, effect.matrix, K, se_m=se, covar = c.mat, bic.var = args$bic_var)
#snp_ids = paste0("rs", 1:nrow(X))
#trait_names = paste0("T",1:ncol(X))
gwasMFGrid <- function(X,W, snp_ids, trait_names, K=0, param.iters = 10, subsample = FALSE)
{
  d <- initializeGwasMF(X,W, snp_ids, trait_names, K=0)
  option <- d$options; args <- d$args; hp <- d$hp; all_ids <- d$all_ids; names <- d$names
  args$auto_grid_search <- TRUE
  option$subsample <- subsample

    #Not sure if I need this now...
  if(FALSE) 
  {
    Sigma <- input.dat$C; 
    #make it a block matrix for this analysis...
    if(is.null(Sigma))
    {
      option$covar <- diag(ncol(X))
    }else
    {
      option$covar <- blockifyCovarianceMatrix(Sigma)
    }
  }
  message(paste0("Number of SNPs for analysis: ", nrow(X)))
  message(paste0("Number of traits for analysis: ", ncol(X)))

  message("initializing for automatic grid search")
  init.params <- GetCoarseSparsityParams(X,W,option, burn.in.reps = 1, logs = TRUE, n.points = 4)
  alphas <- init.params$alphas; lambdas <- init.params$lambdas
  message("Alphas")
  print(alphas)
  message("Lambdas")
  print(lambdas)
  
  grid.solutions <- getGridMatrices(args, option, alphas, lambdas, 7)

  check_stats = max(sapply(grid.solutions, length))
  if(check_stats == 1)
  {
    log_print("No runs completed. Terminating")
    quit()
  }
  
  ret <- GetNextParams(alphas, lambdas,init.k, grid.solutions)
  option$conv0 <- 0.01
  grid.solutions.round.final <- getGridMatrices(args, option, ret$alphas, ret$lambdas, 10)
  
  return(ret)
  }


GetGridPerformanceMetrics <- function(alphas, lambdas,grid.output, init.k)
{
  presumed.order <-unlist(lapply(alphas, function(i) lapply(lambdas, function(j) paste0(i,"_",j))))
  #Fi
  init.k <- 9
  final.out <- NULL
  for(run.i in 1:length(grid.output))
  {
    run <- grid.output[[run.i]]
    #If sparsity is too high or low, 
    v.only <- lapply(run, function(x) x$V)
    u.only <- lapply(run, function(x) x$U)
    if(length(v.only) == 0)
    {
      message("stop here, no returns.")
    }
    coph.results <-  copheneticAssessment(v.only)
    corr.results.per.k2 <- scaledCorrelationAssessment(v.only, type = "squared")
    #The median of the averages
    corr.avg <- scaledCorrelationAssessment(v.only, type = "average_r2")
    nruns <- length(v.only)
    #non-empty entries
    non.empty.count <- nonEmptyAvg(v.only,u.only)
    v.sparsity <- matrixSparsityAvg(v.only,init.k )
    u.sparsity <- matrixSparsityAvg(u.only,init.k)
    v.sparsity.init <- matrixSparsityAvg(v.only, init.k, wrt.init = TRUE)
    u.sparsity.init <- matrixSparsityAvg(u.only, init.k, wrt.init = TRUE)
    #PVE:
    pve <- mean(sapply(run, function(x) sum(x$PVE)))
    
    final.out <- rbind(final.out, 
                       data.frame("Coph" = unlist(coph.results), "Corr_per_factor_pair_SQ"=unlist(corr.results.per.k2), 
                                  "Med_Average_R2" = unlist(corr.avg), 
                                  "Settings" = presumed.order[run.i], "Init_K" = init.k, "Group_nonempty_average" = non.empty.count, 
                                  "V_sparsity_raw" = unlist(v.sparsity), "U_sparsity_raw" = unlist(u.sparsity), 
                                  "V_sparsity_global" = unlist(v.sparsity.init), "U_sparsity_global" = unlist(u.sparsity.init), 
                                  "nruns" = unlist(nruns), "Avg_pve"=unlist(pve)))
    #need the cohenetic correlation across all from the same group
  }
  final.out <- data.frame(final.out %>% tidyr::separate(Settings, into = c("Alpha", "Lambda"), sep= "_",remove = FALSE))
  by.group.med <- final.out %>% group_by(Group_nonempty_average) %>% summarise("Coph_med" = median(Coph))
  left_join(final.out, by.group.med, by = "Group_nonempty_average")
}

  
GetNextParams <- function(alphas, lambdas,init.k, all.solutions)
{
  perf.df <- GetGridPerformanceMetrics(alphas, lambdas,all.solutions, init.k)
  tops <- perf.df %>% tidyr::drop_na() %>% filter(Coph > 0.9) %>% filter(V_sparsity_global > 0.05, V_sparsity_raw > 0.05) %>% 
             arrange(Med_Average_R2) %>% filter(!is.na(Coph_med), Coph_med > 0.9, Avg_pve > 0.5, U_sparsity_raw > 0.05)
    
  #how do we narrow this down to a small number of points to search?
  best.a <- unique(as.numeric(gsub(tops$Alpha, pattern = "A", replacement = "")))
  best.l <- unique(as.numeric(gsub(tops$Lambda, pattern = "L", replacement = "")))
  if(nrow(tops) >= 1)
  {
    new.a <- ProposeSparsityParamsFromGrid(best.a,alphas, 4)
    new.l <- ProposeSparsityParamsFromGrid(c(best.l),lambdas,4)
    print("Alpha:")
    print(paste0((new.a),collapse = ","))
    print("Lambda:")
    print(paste0((new.l),collapse = ","))
  }else #(nrow(tops) == 0)
  {
    message("Parameters too strict. Choosing the most stable solution with any sparsity and K > 1.")
    #Just choosing the most stable one with some sparsity
    tops <- data.frame(final.out) %>% tidyr::drop_na() %>% 
      filter(Group_nonempty_average > 1, V_sparsity_raw > 0) %>% arrange(-Coph) %>% slice(1)  
  }
  return(list("alphas" = new.a, "lambdas" = new.l))
  
}

if(FALSE)
{
  

  option_list <- list(
    make_option(c("--lambda_gc"), type = 'character', help = "Specify the path to the genomic correction coefficients. If none provided, none used", default = ""),
    make_option(c("--trait_names"), type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the orderr in the input tables.", default = ""),
    make_option(c("--calibrate_k"), type = "logical", help = "Give just an estimate of K", action = "store_true", default = FALSE),
    make_option(c("--fixed_first"), type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE),
    make_option(c("--debug"), type = "logical", help = "if want debug run", action = "store_true", default = FALSE),
    make_option(c("--overview_plots"), type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE),
    make_option(c("-k", "--nfactors"), type = "character", help = "specify the number of factors", default = "0"),
    make_option(c("-i", "--niter"), type = "integer", help = "specify the number of iterations", default = 30),
    make_option(c("--posF"), type = "logical", default = FALSE,  help = "Specify if you want to use the smei-nonnegative setup.", action = "store_true"),
    make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand"),
    make_option(c("--epsilon"), type="double", default= 1e-8, help="The convergence criteria for the L1. If exploring the space, try making this larger to speed up runtime. "),
    make_option(c("--subsample"), type="integer", default= 0, help="Specify if you wish to subsample from the full dataset, and if so by how much. For repeated runs, this will choose a different subsample each time."),
    make_option(c("-c","--covar_matrix"), type = "character", default = "", help = "Path to a covariance matrix. Should be pre-filtered to include only significant elements")
  )
}

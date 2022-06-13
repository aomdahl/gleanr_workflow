#Source everything you need:

#... you know, I am pretty sure yuan has a setup for this already. like to specify the desired sparsity or whatever.
pacman::p_load(tidyr, plyr, dplyr, ggplot2, stringr, penalized, cowplot, parallel, doParallel, Xmisc, logr, coop, data.table)
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



parser <- ArgumentParser$new()
parser$add_description("Script to run matrix factorization")
parser$add_argument("--gwas_effects", type = 'character', help = "Specify the Z or B file, depending on specified weighting scheme. First column is ids of each variant, column names specify the trait")
parser$add_argument("--uncertainty", type = 'character', help = "Specify the path to the SE or other uncertainty file, depending on the weightin scheme.irst column is ids of each variant, column names specify the trait")
parser$add_argument("--lambda_gc", type = 'character', help = "Specify the path to the genomic correction coefficients. If none provided, none used", default = "")
parser$add_argument("--trait_names", type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the orderr in the input tables.", default = "")
parser$add_argument("--weighting_scheme", type = 'character', help = "Specify either Z, B, B_SE, B_MAF", default = "B_SE")
parser$add_argument("--init_F", type = 'character', default = "ones_eigenvect", help = "Specify how the F matrix should be initialized. Options are [ones_eigenvect (svd(cor(z))_1), ones_plain (alls 1s), plieotropy (svd(cor(|z|))_1)]")
parser$add_argument("--init_L", type = 'character', help = "Specify this option to start by initializing L rather than F. Options are [random, pleiotropy]. Use empty string for none (default)", default = "")
parser$add_argument("--alphas", type = 'character', help = "Specify which alphas to do, all in quotes, separated by ',' character")
parser$add_argument("--lambdas", type = 'character', help = "Specify which lambdas to do, all in quotes, separated by ',' character")
parser$add_argument("--scaled_sparsity", type = "logical", action = "store_true", default = FALSE, help = "Specify this to scale the sparsity params by dataset to be between 0 and 1")
parser$add_argument("--output", type = "character", help = "Source file location")
parser$add_argument("--cores", type = "numeric", help = "Number of cores", default = 1)
parser$add_argument("--fixed_first", type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE)
parser$add_argument("--debug", type = "logical", help = "if want debug run", action = "store_true", default = FALSE)
parser$add_argument("--overview_plots", type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE)
parser$add_argument("-k", "--nfactors", type = "numeric", help = "specify the number of factors", default = 15)
parser$add_argument("-i", "--niter", type = "numeric", help = "specify the number of iterations", default = 30)
parser$add_argument("--posF", type = "logical", default = FALSE,  help = "Specify if you want to use the smei-nonnegative setup.", action = "store_true")
parser$add_argument("--scale_n", type = "character", default = "",  help = "Specify the path to a matrix of sample sizes if you want to scale by sqrtN as well as W")
parser$add_argument("--IRNT", type = "logical", default = FALSE,  help = "If you want to transform the effect sizes ", action = "store_true")
parser$add_argument("--autofit", type = "character", default = "",  help = "Specify a 1 if you want to autofit sparsity parameter for the whole thing..")
parser$add_argument("-f", "--converged_F_change", type="double", default=0.02,help="Change in factor matrix to call convergence")
parser$add_argument("-o", "--converged_obj_change", type="double", default=1,help="Relative change in the objective function to call convergence")
parser$add_argument("--no_SNP_col", type="logical", default= FALSE, action = "store_true", help="Specify this option if there is NOT first column there...")
parser$add_argument("--parameter_optimization", type="numeric", default= 1, action = "store_true", help="Specify how many iterations to run to get the cophenetic optimization, etc. ")
parser$add_argument("--genomic_correction", type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand")
parser$add_argument("-v", "--verbosity", type="numeric", default= 1, help="How much output information to give in report? 0 is quiet, 1 is loud")

parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
message("Please make sure the first column of input data is SNP/RSIDs.")
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
  
  args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/variant_lists/hm3.pruned.beta.txt"
  args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/variant_lists/hm3.pruned.se.txt"
  args$scale_n <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/variant_lists/hm3.pruned.N.txt"
  args$genomic_correction <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/variant_lists/hm3.pruned.GC.txt"
  args$nfactors <- 6
  args$niter <- 10
  args$alphas <- "0.01,0.05,0.1,0.15,0.2" #using 0.00001 since it doesn't like 0
  args$lambdas <- "0.01,0.05,0.1,0.15,0.2"
args$cores <- 1
args$fixed_first <- TRUE
args$weighting_scheme = "B_SE"
args$output <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/gwasMF_hm3.tmp"
args$converged_obj_change <- 1
args$scaled_sparsity <- TRUE
args$posF <- FALSE
args$autofit <- 0
}

lf <- log_open(paste0(args$output, "gwasMF_log.", Sys.Date(), ".txt"), show_notes = FALSE)

output <- args$output

#Read in the hyperparameters to explore
hp <- readInParamterSpace(args)
log_print("Input sparsity settings")
log_print("Alphas: ")
log_print(hp$a)
log_print("Lambdas: ")
log_print(hp$l)

#set up settings
option <- readInSettings(args)
log_print("Succesfully loaded program options....")
#read in the data
input.dat <- readInData(args)
X <- input.dat$X; W <- input.dat$W; all_ids <- input.dat$ids; names <- input.dat$trait_names

#Estimate sparsity space (doesn't work too well.)
max_sparsity <- approximateSparsity(X, W, option)
log_print("Estimating sparsity maximums based on an SVD approximation.")
log_print(paste0("Max alpha: ", max_sparsity$alpha))
log_print(paste0("Max lambda: ", max_sparsity$lambda))


if(option$calibrate_sparsity)
{
  option$max_alpha <- max_sparsity$alpha
  option$max_lambda <- max_sparsity$lambda
  if(all(hp$a <= 1) & all(hp$a >= 0) & all(hp$l <= 1) & all(hp$l >= 0) & option$calibrate_sparsity)
  {
    alphas <- hp$a * max_sparsity$alpha
    lambdas <- hp$l * max_sparsity$lambda
    
    option$calibrate_sparsity <- FALSE
    updateLog("Scaled amounts are:", option$V)
    updateLog("    Lambdas:", option$V)
    updateLog(paste(round(lambdas, digits = 2)), option$V)
    updateLog("    Alphas:", option$V)
    updateLog(round(alphas, digits = 3), option$V)
  }else
  {
    message('Inputed sparsity parameters must be between 0 and 1 to use the "calibrate_sparsity option".')
    quit()
  }
} else{
  alphas <- hp$a
  lambdas <- hp$l
}


for(opti in 1:args$parameter_optimization) 
{
  message(paste0("On optimization run number: ", opti))
  #Store stats here
  run_stats <- list()
  type_ <- args$weighting_scheme 
  run_stats <- list()
  name_list <- c()
  
  a_plot <- c()
  l_plot <- c()
  iter_count <- 1
  #In the case of autofit, you wouldn't go through a bunch of starting points, would you? maybe. Actually not a bad idea, huh.
  
  for(i in 1:length(alphas)){
    for (j in 1:length(lambdas)){
      a <- alphas[i]
      l <- lambdas[j]
      option[['alpha1']] <- as.numeric(a)
      option[['lambda1']] <- as.numeric(l)
      option[["parallel"]] <- FALSE
      start <- Sys.time()
      run <- Update_FL(as.matrix(X), as.matrix(W), option)
      time <-  round(difftime(Sys.time(), start, units = "mins"), digits = 3)
      run_stats[[iter_count]] <- c(run, time)
      
      iter_count <- iter_count + 1
      #things for plots
      if(length(run) == 0)
      {
        fname = paste0(output, "A", a, "_L", l, "_","K", args$ type_, ".NO_OUTPUT.png")
      } else
      {
        updateStatement(l,a,hp$l[j],hp$a[i], run,time)
        fname = paste0(output, "A", round(a, digits = 4), "_L", round(l, digits = 4), "_", type_,".", opti, ".png")
        title_ = paste0("A", round(a, digits = 4), "_L", round(l, digits = 4), " Type = ", args$weighting_scheme )
        p <- plotFactors(run[[1]],trait_names = names, title = title_, scale.cols = TRUE)
        ggsave(filename = fname, plot = p, device = "png", width = 10, height = 8)
        #Lazy bookeeping nonsense for plots
        name_list <- c(name_list,  paste0("A", round(a, digits = 6), "_L", round(l, digits = 6)))
        a_plot <- c(a_plot, a)
        l_plot <- c(l_plot, l)
        writeFactorMatrix(names, all_ids, run,  paste0("A", round(a, digits = 6), "_L", round(l, digits = 6), "_", type_, ".", opti), output)
      } 
      
    }
  }
  #Save all the data...
  #We actually need output information
  save(run_stats, file = paste0(output, "runDat.", opti, ".RData"))
  log_print(paste0("Run data written out to: ", output, "runDat.RData"))
  #writeFactorMatrices(round(alphas, digits = 3), round(lambdas, digits = 3), names,all_ids, run_stats,output)
  check_stats = max(sapply(run_stats, length))
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
    
   # objectives <- data.frame(lapply(run_stats, function(x) x[[6]])) %>% drop_na() 
   # colnames(objectives) <- name_list
   # objectives <- objectives %>% mutate("iteration" =1:nrow(objectives)) %>% reshape2::melt(., id.vars = "iteration")
   # ggplot(data = objectives, aes(x = iteration, y = value, color = variable)) + geom_point() + geom_line() + 
   #   theme_minimal(15) + ggtitle("Objective function") + ylab("Objective score")
   # ggsave(paste0(output, "objective.", opti, ".png"), device = "png")    
  }
  rm(run_stats)
  log_print(paste0("Completed iteration number ", opti))
  log_print("----------------------------------")
  log_print("")
}

####Propose the optimal solution
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

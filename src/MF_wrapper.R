#Source everything you need:

#... you know, I am pretty sure yuan has a setup for this already. like to specify the desired sparsity or whatever.
pacman::p_load(tidyr, dplyr, ggplot2, stringr, penalized, cowplot, parallel, doParallel, Xmisc, logr)
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_F.R")
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/update_FL.R")
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_L.R")
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/plot_functions.R")
source('/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/compute_obj.R')
source('/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/buildFactorMatrices.R')
source('/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/sparsity_scaler.R')

quickSort <- function(tab, col = 1)
{
  tab[order(tab[,..col], decreasing = TRUE),]
}

updateStatement  <- function(l,a,l_og, a_og, run,time)
{
  
  log_print(paste0("Run complete. (Time: ", time, ")"))
  log_print(paste0("Sparsity params: F ", l, ", L ", a))  
  log_print(paste0("Sparsity scale: F ", l_og, ", L ", a_og))  
  log_print(paste0("Current F sparsity: ", run$F_sparsity))  
  log_print(paste0("Current L sparsity: ", run$L_sparsity))  
  log_print(paste0("Number of active factors:", ncol(run$F)))
}


parser <- ArgumentParser$new()
parser$add_description("Script to run matrix factorization")
parser$add_argument("--gwas_effects", type = 'character', help = "Specify the Z or B file, depending on specified weighting scheme. First column is ids of each variant, column names specify the trait")
parser$add_argument("--uncertainty", type = 'character', help = "Specify the path to the SE or other uncertainty file, depending on the weightin scheme.irst column is ids of each variant, column names specify the trait")
parser$add_argument("--trait_names", type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the orderr in the input tables.", default = "")
parser$add_argument("--weighting_scheme", type = 'character', help = "Specify either Z, B, B_SE, B_MAF", default = "B_SE")
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
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
message("Please make sure the first column of input data is SNP/RSIDs.")
lf <- log_open(paste0(args$output, "/gwasMF_log.", Sys.Date(), ".txt"))
#TODO:
    #add in functionality in the event the objective begins to increase again. Don't want that....
  #finesse functionality to allow for trait-specific variance to be learned (?)
#easy enough to add in.
if(FALSE) #For debug functionality
{
  args <- list()
  args$gwas_effects <- "/work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/beta_signed_matrix.tsv"
  args$uncertainty <- "/work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/se_matrix.tsv"
  args$nfactors <- 6
  args$niter <- 10
  args$alphas <- "1e-16,0.001,0.005,0.01" #using 0.00001 since it doesn't like 0
  args$lambdas <- "1e-16,0.001,0.005,0.01"

args$cores <- 1
args$fixed_first <- TRUE
args$weighting_scheme = "B_SE"
args$output <- "/work-zfs/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/gwasMF_"
args$converged_obj_change <- 1
args$scaled_sparsity <- TRUE
args$posF <- FALSE
args$scaled_sparsity <- TRUE
}
lf <- log_open(paste0(args$output, "/gwasMF_log.", Sys.Date(), ".txt"))


#Read in the hyperparameters to explore
alphas_og <- as.numeric(scan(text = args$alphas, what = character(), sep = ',', quiet = TRUE))
lambdas_og <- as.numeric(scan(text = args$lambdas, what = character(), sep = ',', quiet = TRUE))
output <- args$output
log_print("Input sparsity settings")
log_print(paste0("Alphas: ", alphas))
log_print(paste0("Lambdas: ", lambdas))


#Load the effect size data
effects <- fread(args$gwas_effects) %>% drop_na()
effects <- effects[order(effects[,1], decreasing = TRUE),]
all_ids <- unlist(effects[,1])


#Get the trait names out
if(args$trait_names == "")
{
  message("No trait names provided. Using the identifiers in the tabular effect data instead.")
  names <- names(effects)[-1]
} else{
    names <- scan(args$trait_names, what = character(), quiet = TRUE)
}

effects <- effects[,-1]
if(args$IRNT)
{
  library(RNOmni)
  effects <- apply(effects, 2, RankNorm)
}



#Look at the weighting scheme options...
if(args$weighting_scheme == "Z" || args$weighting_scheme == "B")
{
  message("No scaling by standard error will take place. Input to uncertainty being ignored.")
  W <- matrix(1,nrow = nrow(z), ncol = ncol(z))
  X <- effects
  
} else if(args$weighting_scheme == "B_SE")
{
  message("Scaling by 1/SE.")
  W_se <- fread(args$uncertainty) %>% drop_na() %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.)
  stopifnot(!any(all_ids != W_se[,1]))
  W_se <- W_se %>% select(-1)
  W <- 1/ W_se
  X <- effects
  
} else if(args$weighting_scheme == "B_MAF")
{
  message("Scaling by 1/var(MAF)")
  W_maf <- fread(args$uncertainty) %>% drop_na() %>% 
    filter(ids %in% all_ids) %>% arrange(ids) %>% select(-ids) 
  W <- 1/matrix(apply(W_maf, 2, function(x) 2*x*(1-x)), nrow = nrow(W_maf), ncol = ncol(W_maf))
  X <- effects 
} else
{
  message("No form selected. Please try again.")
  quit()
}

if(args$scale_n != "")
{
  N <- fread(args$scale_n) %>% drop_na() %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.)
  stopifnot(!any(all_ids != N[,1]))
  N <- N %>% select(-1)
  W <- W * (1/sqrt(N))
}

#set up settings
option <- list()
option[['K']] <- args$nfactors
option[['iter']] <- args$niter
option[['convF']] <- 0
option[['convO']] <- args$converged_obj_change
option[['ones_mixed_std']] <- FALSE
option[['ones_mixed_ref']] <- FALSE
option[['ones']] <- FALSE
option[['disp']] <- FALSE
option[['convF']] <- 0
option[['convO']] <- 0
option[['ones_mixed_std']] <- FALSE
option[["preinitialize"]] <- FALSE
option[["ones_mixed_ref"]] <- FALSE
option[['ones_eigenvect']] <- TRUE
option[['ones_plain']] <- FALSE
option[['reweighted']] <- FALSE
option[["glmnet"]] <- FALSE
option[["parallel"]] <- FALSE
option[["fastReg"]] <- FALSE
option[["ridge_L"]] <- FALSE
option[["posF"]] <- args$posF
option[["autofit"]] <- args$autofit
option$intercept_ubiq <- FALSE
option$traitSpecificVar <- FALSE
option$calibrate_sparsity <- args$scaled_sparsity
if(args$cores > 1)
{
  message(paste("Running in parallel on", args$cores, "cores"))
  option[["parallel"]] <- TRUE
}
option[["ncores"]] <- args$cores
option[["fixed_ubiq"]] <- args$fixed_first
#Store stats here
run_stats <- list()
type_ <- args$weighting_scheme 
run_stats <- list()
name_list <- c()

a_plot <- c()
l_plot <- c()
iter_count <- 1


#Actually run the factorization

#recommend
#3/28 attempt

#max_sparsity <- Update_FL(as.matrix(X), as.matrix(W), option)
max_sparsity <- approximateSparsity(as.matrix(X), as.matrix(W), option)

log_print("Scaling sparsity by an SVD-based max.")
log_print(paste0("Max alpha: ", max_sparsity$alpha))
log_print(paste0("Max lambda: ", max_sparsity$lambda))


if(option$calibrate_sparsity)
{
  option$max_alpha <- max_sparsity$alpha
  option$max_lambda <- max_sparsity$lambda
  if(all(alphas_og <= 1 & alphas_og >= 0 & lambdas_og <= 1 & lambdas_og >= 0) & option$calibrate_sparsity)
  {
    alphas <- alphas_og * max_sparsity$alpha
    
    lambdas <- lambdas_og * max_sparsity$lambda
    
    option$calibrate_sparsity <- FALSE
    message("Scaled amounts are:")
    message("    Lambdas:")
    message(round(lambdas, digits = 3))
    message("    Alphas:")
    message(round(alphas, digits = 3))
  }else
  {
    message('Inputed sparsity parameters must be between 0 and 1 to use the "calibrate_sparsity option".')
    quit()
  }
} else{
  alphas <- alphas_og
  lambdas <- lambdas_og
}
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
    end <- Sys.time()
    time <- end-start
    run_stats[[iter_count]] <- c(run, time)
    
    iter_count <- iter_count + 1
    #things for plots
    if(length(run) == 0)
    {
        fname = paste0(output, "A", a, "_L", l, "_","K", args$ type_, ".NO_OUTPUT.png")
    } else
    {
      updateStatement(l,a,lambdas_og[j],alphas_og[i], args$run,time)
      fname = paste0(output, "A", round(a, digits = 3), "_L", round(l, digits = 3), "_", type_, ".png")
        title_ = paste0("A", round(a, digits = 3), "_L", round(l, digits = 3), " Type = ", args$weighting_scheme )
          p <- plotFactors(run[[1]],trait_names = names, title = title_)
          ggsave(filename = fname, plot = p, device = "png", width = 10, height = 8)
          #Lazy bookeeping nonsense for plots
          name_list <- c(name_list,  paste0("A", round(a, digits = 3), "_L", round(l, digits = 3)))
          a_plot <- c(a_plot, a)
          l_plot <- c(l_plot, l)   
    } 

  }
}
  #Save all the data...
#We actually need output information
save(run_stats, file = paste0(output, "runDat.RData"))
log_print(paste0("Run data written out to: ", output, "runDat.RData"))
writeFactorMatrices(round(alphas, digits = 3), round(lambdas, digits = 3), names,all_ids, run_stats,output)
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
  message(paste0(output, "factor_sparsity.png"))

  ggplot(data = factor_sparsities, aes(x = alpha, y= lambda, fill = sparsity)) + geom_tile() + 
    theme_minimal(15) + scale_fill_gradient(low="white", high="navy") + ggtitle("Factor Sparsity")
  ggsave(paste0(output, "factor_sparsity.png"), device = "png")  

  ggplot(data = loading_sparsities, aes(x = alpha, y= lambda, fill = sparsity)) + geom_tile() + 
    theme_minimal(15) + scale_fill_gradient2(low="white", high="navy", limits = c(0,max(loading_sparsities$sparsity))) + ggtitle("Loading Sparsity")
  ggsave(paste0(output, "loading_sparsity.png"), device = "png")  

  #Plot the change in objective function for each one too....

  objectives <- data.frame(lapply(run_stats, function(x) x[[6]])) %>% drop_na() 
  colnames(objectives) <- name_list
  objectives <- objectives %>% mutate("iteration" =1:nrow(objectives)) %>% reshape2::melt(., id.vars = "iteration")
  ggplot(data = objectives, aes(x = iteration, y = value, color = variable)) + geom_point() + geom_line() + 
    theme_minimal(15) + ggtitle("Objective function") + ylab("Objective score")
  ggsave(paste0(output, "objective.png"), device = "png")    
}
log_print("Option settings",console = FALSE)
log_print(option,console = FALSE)
log_close()
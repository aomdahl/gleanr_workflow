#Source everything you need:

#... you know, I am pretty sure yuan has a setup for this already. like to specify the desired sparsity or whatever.
pacman::p_load(tidyr, dplyr, readr, ggplot2, stringr, penalized, cowplot, parallel, doParallel, Xmisc)
source("src/fit_F.R")
source("src/update_FL.R")
source("src/fit_L.R")
source("src/plot_functions.R")
source('src/compute_obj.R')
source('src/buildFactorMatrices.R')

parser <- ArgumentParser$new()
parser$add_description("Script to run matrix factorization")
parser$add_argument("--gwas_effects", type = 'character', help = "Specify the Z or B file, depending on specified weighting scheme. First column is ids of each variant, column names specify the trait")
parser$add_argument("--uncertainty", type = 'character', help = "Specify the path to the SE or other uncertainty file, depending on the weightin scheme.irst column is ids of each variant, column names specify the trait")
parser$add_argument("--trait_names", type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the orderr in the input tables.")
parser$add_argument("--weighting_scheme", type = 'character', help = "Specify either Z, B, B_SE, B_MAF", default = "B_SE")
parser$add_argument("--alphas", type = 'character', help = "Specify which alphas to do, all in quotes, separated by ',' character")
parser$add_argument("--lambdas", type = 'character', help = "Specify which lambdas to do, all in quotes, separated by ',' character")
parser$add_argument("--output", type = "character", help = "Source file location")
parser$add_argument("--cores", type = "numeric", help = "Number of cores", default = 2)
parser$add_argument("--fixed_first", type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE)
parser$add_argument("--debug", type = "logical", help = "if want debug run", action = "store_true", default = FALSE)
parser$add_argument("--overview_plots", type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE)
parser$add_argument("-k", "--nfactors", type = "numeric", help = "specify the number of factors", default = 15)
parser$add_argument("-i", "--niter", type = "numeric", help = "specify the number of iterations", default = 30)
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()

#TODO:
    #add in functionality in the event the objective begins to increase again. Don't want that....
#Read in the hyperparameters to explore
alphas <- scan(text = args$alphas, what = character(), sep = ',', quiet = TRUE)
lambdas <- scan(text = args$lambdas, what = character(), sep = ',', quiet = TRUE)
output <- args$output
#Load the data
effects <- fread(args$gwas_effects) %>% drop_na() %>% arrange(ids)
all_ids <- effects$ids
if(args$trait_names == "")
{
  message("No trait names provided. Using the identifiers in the tabular effect data instead.")
  names <- names(effects)[-1]
} else{
    names <- scan(args$trait_names, what = character(), quiet = TRUE)
}
effects <- effects %>% select(-ids)
if(args$weighting_scheme == "Z" || args$weighting_scheme == "B")
{
  message("No scaling by standard error will take place. Input to uncertainty being ignored.")
  W <- matrix(1,nrow = nrow(z), ncol = ncol(z))
  X <- effects
  
} else if(args$weighting_scheme == "B_SE")
{
  message("Scaling by 1/SE.")
  W_se <- fread(args$uncertainty) %>% drop_na() %>% filter(ids %in% all_ids) %>% arrange(ids) %>% select(-ids)
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

#set up settings
option <- list()
option[['K']] <- args$nfactors
option[['iter']] <- args$niter
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
if(args$debug)
{
  message(alphas)
  message(lambdas)
  message(head(W))
  message(head(X))
  message("Debug okay?")
  quit()
}
a_plot <- c()
l_plot <- c()
iter_count <- 1
#Actually run the factorization
for(a in alphas){
  for (l in lambdas){
    option[['alpha1']] <- as.numeric(a)
    option[['lambda1']] <- as.numeric(l)
    option[["parallel"]] <- TRUE
    start <- Sys.time()
    run <- Update_FL(X, W, option)
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
        fname = paste0(output, "A", a, "_L", l, "_", type_, ".png")   
    } 
    title_ = paste0("A", a, "_L", l, " Type = ", args$weighting_scheme )
    p <- plotFactors(run[[1]],trait_names = names, title = title_)
    ggsave(filename = fname, plot = p, device = "png", width = 10, height = 8)
    #Lazy bookeeping nonsense for plots
    name_list <- c(name_list,  paste0("A", a, "_L", l))
    a_plot <- c(a_plot, a)
    l_plot <- c(l_plot, l)
  }
}
#Save all the data...
#We actually need output information
save(run_stats, file = paste0(output, "/runDat.RData"))
writeFactorMatrices(alphas, lambdas, names, run_stats,output)

if(args$overview_plots)
{
  factor_sparsities <- data.frame("alpha" = a_plot, "lambda" = l_plot, 
                                  "sparsity" = unlist(lapply(run_stats, function(x) x[[4]])), 
                                  "runtime" = unlist(lapply(run_stats, function(x) x[[7]])))
  loading_sparsities <- data.frame("alpha" = a_plot, "lambda" = l_plot, 
                                  "sparsity" = unlist(lapply(run_stats, function(x) x[[3]])), 
                                  "runtime" = unlist(lapply(run_stats, function(x) x[[7]])))
  message(paste0(output, "/factor_sparsity.png"))

  ggplot(data = factor_sparsities, aes(x = alpha, y= lambda, fill = sparsity)) + geom_tile() + 
    theme_minimal(15) + scale_fill_gradient(low="white", high="navy") + ggtitle("Factor Sparsity")
  ggsave(paste0(output, "/factor_sparsity.png"), device = "png")  

  ggplot(data = loading_sparsities, aes(x = alpha, y= lambda, fill = sparsity)) + geom_tile() + 
    theme_minimal(15) + scale_fill_gradient2(low="white", high="navy", limits = c(0,max(loading_sparsities$sparsity))) + ggtitle("Loading Sparsity")
  ggsave(paste0(output, "/loading_sparsity.png"), device = "png")  
  #Make a plot with the various runtimes.
  ggplot(data = factor_sparsities, aes(x = alpha, y= lambda, fill = runtime)) + geom_tile() + 
    theme_minimal(15) + scale_fill_gradient(low="white", high="navy") + ggtitle("Runtime")
  ggsave(paste0(output, "/runtimes.png"), device = "png")      
  #Plot the change in objective function for each one too....

  objectives <- data.frame(lapply(run_stats, function(x) x[[6]])) %>% drop_na() 
  colnames(objectives) <- name_list
  objectives <- objectives %>% mutate("iteration" =1:nrow(objectives)) %>% reshape2::melt(., id.vars = "iteration")
  ggplot(data = objectives, aes(x = iteration, y = value, color = variable)) + geom_point() + geom_line() + 
    theme_minimal(15) + ggtitle("Objective function") + ylab("Objective score")
  ggsave(paste0(output, "/objective.png"), device = "png")    
}






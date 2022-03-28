#Source everything you need:

#... you know, I am pretty sure yuan has a setup for this already. like to specify the desired sparsity or whatever.
pacman::p_load(tidyr, dplyr, ggplot2, stringr, penalized, cowplot, parallel, doParallel, Xmisc)
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_F.R")
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/update_FL.R")
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/fit_L.R")
source("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/plot_functions.R")
source('/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/compute_obj.R')
source('/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/src/buildFactorMatrices.R')


quickSort <- function(tab, col = 1)
{
  tab[order(tab[,..col], decreasing = TRUE),]
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
parser$add_argument("--cores", type = "numeric", help = "Number of cores", default = 2)
parser$add_argument("--fixed_first", type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE)
parser$add_argument("--debug", type = "logical", help = "if want debug run", action = "store_true", default = FALSE)
parser$add_argument("--overview_plots", type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE)
parser$add_argument("-k", "--nfactors", type = "numeric", help = "specify the number of factors", default = 15)
parser$add_argument("-i", "--niter", type = "numeric", help = "specify the number of iterations", default = 30)
parser$add_argument("--posF", type = "logical", default = FALSE,  help = "Specify if you want to use the smei-nonnegative setup.", action = "store_true")
parser$add_argument("--scale_n", type = "character", default = "",  help = "Specify the path to a matrix of sample sizes if you want to scale by sqrtN as well as W")
parser$add_argument("-f", "--converged_F_change", type="double", default=0.02,help="Change in factor matrix to call convergence")
parser$add_argument("-o", "--converged_obj_change", type="double", default=1,help="Relative change in the objective function to call convergence")
parser$add_argument("--no_SNP_col", type="logical", default= FALSE, action = "store_true", help="Specify this option if there is NOT first column there...")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
message("Please make sure the first column of input data is SNP/RSIDs.")
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
#Read in the hyperparameters to explore
alphas <- as.numeric(scan(text = args$alphas, what = character(), sep = ',', quiet = TRUE))
lambdas <- as.numeric(scan(text = args$lambdas, what = character(), sep = ',', quiet = TRUE))
output <- args$output

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




#Look at the weighting scheme options...
if(args$weighting_scheme == "Z" || args$weighting_scheme == "B")
{
  message("No scaling by standard error will take place. Input to uncertainty being ignored.")
  W <- matrix(1,nrow = nrow(z), ncol = ncol(z))
  X <- effects
  
} else if(args$weighting_scheme == "B_SE")
{
  message("Scaling by 1/SE.")
  W_se <- fread(args$uncertainty) %>% drop_na() %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.) %>% select(-1)
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
  N <- fread(args$scale_n) %>% drop_na() %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.) %>% select(-1)
  W <- W * 1/sqrt(N)
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

#3/24
#If you want the sparsity to be between 0 and 1 always, trying this
if(args$scaled_sparsity & FALSE)
{
  Y_ <- as.matrix(X) * as.matrix(W)
  pca <- svd(Y_)
  lambda_max_L <- t(pca$u) %*% Y_ #equiviland to diag(d) %*% t(V)
  reg_norm_f <- norm(lambda_max_L, "I")#This goes to lambda
  
  lambda_max_f <- t(diag(pca$d) %*% pca$v) %*% t(Y_) 
  reg_norm_l <- norm(lambda_max_f, "I") #this goes to alpha
  
  #update the alphas and others..
  alphas <- as.numeric(alphas) * reg_norm_l
  lambdas <- as.numeric(lambdas) * reg_norm_f
  #tests:
  #Sparsity of 1 should really kill everything. If there is leftover, its not high enough on either side...
  #SParsity of 0 should have none
  #this didn't work- we didn't quite achieve the max, and we were definately too sparse....
  #I could try with +1, and then run. It should short out if its too big though!
  max( abs(t(Y - mean(Y)*(1-mean(Y))) %*% X ) )/ ( alpha * n) # largest lambda value
  max(abs(((Y_ - mean(Y_)*(1-mean(Y_))) %*% lambda_max_L) / dim(Y_)[1]))
  out <-c()
  for(i in 1:93)
  {
    out <- c(out, glmnet(x=lambda_max_L,y=Y_[i,],alpha = 1,standardize=FALSE)$lambda[1])
  }
  
  #not the case my frined. 11.25319 is too small. Way too small.
}


#Actually run the factorization

#recommend
#3/28 attempt

max_sparsity <- Update_FL(as.matrix(X), as.matrix(W), option)
if(option$calibrate_sparsity)
{
  if(all(alphas <= 1 & alphas >= 0 & lambdas <= 1 & lambdas >= 0) & option$calibrate_sparsity)
  {
    alphas <- alphas * max_sparsity$alpha
    
    lambdas <- lambdas * max_sparsity$lambda
    
    option$calibrate_sparsity <- FALSE
    message("Scaled amounts are:")
    message("    Lambdas:")
    message(lambdas)
    message("    Alphas:")
    message(alphas)
  }else
  {
    message('Inputed sparsity parameters must be between 0 and 1 to use the "calibrate_sparsity option".')
    quit()
  }
}

for(a in alphas){
  for (l in lambdas){
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
      message(paste0("Sparsity params: F ", l, ", L ", a))  
      message(paste0("Current F sparsity: ", run$F_sparsity))  
      message(paste0("Current L sparsity: ", run$L_sparsity))  
      message(paste0("Number of active factors:", ncol(run$F)))
      fname = paste0(output, "A", round(a, digits = 3), "_L", round(l, digits = 3), "_", type_, ".png")
        title_ = paste0("A", round(a, digits = 3), "_L", round(l, digits = 3), " Type = ", args$weighting_scheme )
          p <- plotFactors(run[[1]],trait_names = names, title = title_)
          print(plotFactors(run[[1]],trait_names = names, title = title_))
          ggsave(filename = fname, plot = p, device = "png", width = 10, height = 8)
          #Lazy bookeeping nonsense for plots
          name_list <- c(name_list,  paste0("A", a, "_L", l))
          a_plot <- c(a_plot, a)
          l_plot <- c(l_plot, l)   
    } 

  }
}
  #Save all the data...
#We actually need output information
save(run_stats, file = paste0(output, "runDat.RData"))
writeFactorMatrices(alphas, lambdas, names,all_ids, run_stats,output)
check_stats = max(sapply(run_stats, length))
print(check_stats)
print(run_stats)
if(check_stats == 1)
{
  message("No runs completed. Terminating")
  quit()
}
if(args$overview_plots)
{
  factor_sparsities <- data.frame("alpha" = a_plot, "lambda" = l_plot, 
                                  "sparsity" = unlist(lapply(run_stats, function(x) x[[4]])), 
                                  "runtime" = unlist(lapply(run_stats, function(x) x[[7]])))
  loading_sparsities <- data.frame("alpha" = a_plot, "lambda" = l_plot, 
                                  "sparsity" = unlist(lapply(run_stats, function(x) x[[3]])), 
                                  "runtime" = unlist(lapply(run_stats, function(x) x[[7]])))
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






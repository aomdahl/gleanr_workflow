#Source everything you need:

#... you know, I am pretty sure yuan has a setup for this already. like to specify the desired sparsity or whatever.
pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, penalized, cowplot, parallel, doParallel, Xmisc)
source("src/fit_F.R")
source("src/update_FL.R")
source("~src/fit_L.R")
source("src/plot_functions.R")
source('src/compute_obj.R')

parser <- ArgumentParser$new()
parser$add_description("Script to explore different weighting schemes")
parser$add_argument("--weighting_scheme", type = 'character', help = "Specify either Z, B, B_SE, B_MAF")
parser$add_argument("--alphas", type = 'character', help = "Specify which alphas to do, all in quotes, separated by ',' character")
parser$add_argument("--lambdas", type = 'character', help = "Specify which lambdas to do, all in quotes, separated by ',' character")
parser$add_argument("--output", type = "character", help = "Source file location")
parser$add_argument("--cores", type = "numeric", help = "Number of cores", default = 2)
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
if(FALSE)
{
  args <- list()
  args$alphas <- "1,2,3"
  args$lambdas <- "1,2,3"
  args$weighting_scheme <- "B_SE"
}
alphas <- scan(text = args$alphas, what = character(), sep = ',')
lambdas <- scan(text = args$lambdas, what = character(), sep = ',')

#Load the data
dir <- "/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/"
z_scores <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.z.tsv") %>% drop_na() %>% arrange(ids)
names <- scan("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.names.tsv", what = character())
z <- as.matrix(z_scores %>% select(-ids))

betas <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.beta.tsv") %>% 
  drop_na() %>% filter(ids %in% z_scores$ids) %>% arrange(ids) %>% select(-ids) %>% as.matrix() 
all_ids <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.beta.tsv") %>% 
  drop_na() %>% filter(ids %in% z_scores$ids) %>% arrange(ids) %>% select(ids)

rsid_map <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5.pruned_rsids.txt", header =FALSE) %>% rename("ids" = V1, "rsids" = V2)
if(args$weighting_scheme == "Z")
{
  W <- matrix(1,nrow = nrow(z), ncol = ncol(z))
  X <- z
  
} else if(args$weighting_scheme == "B"){
  W <- matrix(1,nrow = nrow(z), ncol = ncol(z))
  X <- betas
  
} else if(args$weighting_scheme == "B_SE")
{
  W_se <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.se.tsv") %>% drop_na() %>% filter(ids %in% z_scores$ids) %>% arrange(ids) %>% select(-ids) %>% as.matrix()
  W <- 1/ W_se
  X <- betas
  
} else if(args$weighting_scheme == "B_MAF")
{
  W_maf <- fread(paste0(dir, "seed2_thresh0.9_h2-0.1_vars1e-5.FULL_LIST.1e-5.maf.tsv") %>% 
    drop_na() %>% filter(ids %in% z_scores$ids) %>% arrange(ids) %>% select(-ids) %>% as.matrix()
  W <- 1/matrix(apply(W_maf, 2, function(x) 2*x*(1-x)), nrow = nrow(W_maf), ncol = ncol(W_maf))
  X <- betas
} else
{
  print("No form selected. Please try again.")
  quit()
}



#set up settings
option <- list()
option[['K']] <- 15
option[['alpha1']] <- 1
option[['lambda1']] <- 1
option[['iter']] <- 30
option[['ones_mixed_std']] <- FALSE
option[['ones_mixed_ref']] <- FALSE
option[['ones']] <- FALSE
option[["fixed_ubiq"]] <- FALSE
option[['disp']] <- FALSE
option[['convF']] <- 0
option[['convO']] <- 0
option[['ones_mixed_std']] <- FALSE
option[["ones_mixed_ref"]] <- FALSE
option[['ones_eigenvect']] <- TRUE
option[['ones_plain']] <- FALSE
option[['reweighted']] <- FALSE
option[["glmnet"]] <- FALSE
option[["parallel"]] <- TRUE
option[["ncores"]] <- args$cores
#Store stats here
run_stats <- list()


f_sparsity <- 0
target <- 0.87
iter_count = 1
type = args$weighting_scheme 
run_stats <- list()
for(a in alphas){
  for (l in lambdas){
    option[['alpha1']] <- a
    option[['lambda1']] <- l
    option[["parallel"]] <- TRUE
    start <- Sys.time()
    run <- Update_FL(X, W, option)
    print(a)
    print(l)
    end <- Sys.time()
    time <- end-start
    run_stats[[iter_count]] <- c(run, time)
    iter_count <- iter_count + 1
    title_ = paste0("A", a, "_L", l, "Type = ", args$weighting_scheme )
    print(title_)
    p <- plotFactors(run[[1]],trait_names = names, title = title_)
    ggsave(paste0(p, output, "_A", a, "_L", l, "_", type, ".png"))
  }
}

#Make a plot with the various sparsities
#alpha <- 1:9
#lambda <- 1:9
factor_sparsities <- data.frame("alpha" = alpha, "lambda" = lambda, 
                                "sparsity" = unlist(lapply(run_stats, function(x) x[[4]])), 
                                "runtime" = unlist(lapply(run_stats, function(x) x[[7]])))
loading_sparsities <- data.frame("alpha" = alpha, "lambda" = lambda, 
                                 "sparsity" = unlist(lapply(run_stats, function(x) x[[3]])), 
                                 "runtime" = unlist(lapply(run_stats, function(x) x[[7]])))


ggplot(data = factor_sparsities, aes(x = alpha, y= lambda, color = sparsity)) + geom_point() + 
  theme_minimal(15) + scale_color_gradient(low="navy", high="red") + ggtitle("Factor Sparsity")
ggsave(paste0(output, "/factor_sparsity.png"))                                                         
ggplot(data = loading_sparsities, aes(x = alpha, y= lambda, color = sparsity)) + geom_point() + 
  theme_minimal(15) + scale_color_gradient(low="navy", high="red") + ggtitle("Loading Sparsity")
ggsave(paste0(output, "/loading_sparsity.png"))  
#Make a plot with the various runtimes.
ggplot(data = factor_sparsities, aes(x = alpha, y= lambda, color = runtime)) + geom_point() + 
  theme_minimal(15) + scale_color_gradient(low="navy", high="red") + ggtitle("Runtime")
ggsave(paste0(output, "/runtimes.png"))  
#Plot the change in objective function for each one too....
#run_stat <- run_stats_l1
objectives <- data.frame(lapply(run_stats, function(x) x[[6]])) %>% drop_na() %>% mutate("iteration" =1: (option[['iter']]+1))%>%
  reshape2::melt(., id.vars = "iteration")
ggplot(data = objectives, aes(x = iteration, y = value, color = variable)) + geom_point() + geom_line() + 
  theme_minimal(15) + ggtitle("Objective function") + ylab("Objective score")
ggsave(paste0(output, "/objective.png"))  
#Save all the data...
save(run_stats, file = paste0(output, "/runDat.RData"))





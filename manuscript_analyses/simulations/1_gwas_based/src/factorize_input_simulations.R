library(optparse)
#pacman::p_load(data.table, doParallel, foreach, tidyr, dplyr)
pacman::p_load(data.table, tidyr, dplyr)



stdPCAEvals <- function(sim.zscores,prefix="")
{
  ret.list <- list()
  message("Performing standard SVD on all of them:")
  ret.list[[paste0(prefix, "sim.pca.std")]] <- pcaOnAll(sim.zscores, cols = "ALL")
  
  message("Performing BiocSinglar PCA for model selection tools.")
  ret.list[[paste0(prefix, "sim.pca.alt")]] <- lapply(sim.zscores, function(x) BiocSingular::runPCA(x=x, rank=rank-1, center=FALSE))
  
  #for the eblow: this above gives the sds
  #proportion var based on the PCAtools package documentation
  message("Performing model selection using basic heuristics")
  total.vars <- lapply(sim.zscores, function(x) sum(matrixStats::colVars(x)))
  ret.list[[paste0(prefix, "proportionvar")]] <- lapply(1:length(ret.list[[paste0(prefix, "sim.pca.alt")]]), function(i) ret.list[[paste0(prefix, "sim.pca.alt")]][[i]]$sdev^2/total.vars[[i]]*100) #from the PCA code
  ret.list[[paste0(prefix, "gd.point")]] <- lapply(1:length(ret.list[[paste0(prefix, "proportionvar")]]), 
                               function(i) PCAtools::chooseGavishDonoho(sim.zscores[[i]], 
                                                                        var.explained=ret.list[[paste0(prefix, "proportionvar")]][[i]],
                                                                        noise=1))
  ret.list[[paste0(prefix, "elbow.point")]] <- lapply(ret.list[[paste0(prefix, "proportionvar")]], function(x) PCAtools::findElbowPoint(x))
  ret.list[[paste0(prefix, "kaiser.point")]] <- sapply(ret.list[[paste0(prefix, "sim.pca.std")]], function(x) sum(x$d^2 > mean(x$d^2))) #not actually kaiser
  ret.list
}










option_list <- list(
  make_option(c("--input_dir"), type = 'character', help = "Specifgy the input directory with the files- assuming the newest version"),
  make_option(c("--with_covar"), type = 'logical', help = "Use flag if want to include covariance adjustment. False by default", default = FALSE, action="store_true"),
  make_option(c("-o", "--output"), type = 'character', help = "Output file path"),
  make_option(c("-r", "--run_only"), type = 'character', help = "Specify a list of specific methods you want to run. added for later utility.", default = "ALL"),
  make_option(c("-s", "--start"), type = 'integer', help = "Start index", default = 1),
  make_option(c("-c", "--ncores"), type = 'numeric', help = "how many cores you want. if not specified, detcts (dangerous)", default=-1),
  make_option(c("--convergence_criteria"), type = 'numeric', help = "Convergence stopping setting for gleaner", default = 0.001)
)

quiet <- function(x) { 
  sink(sink('/dev/null')) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

odir="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/N30000_multi-overlap_SD-1_h2-0-0_bmi-like.SCALED_Aug2024/"

t <- c("--input_dir=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/N30000_multi-overlap_SD-1_h2-0-0_bmi-like.SCALED_Aug2024//",
       "--with_covar", "--convergence_criteria=0.01",
       "--output=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/N30000_multi-overlap_SD-1_h2-0-0_bmi-like.SCALED_Aug2024//factorization.summaries.RData")
#argv <- parse_args(OptionParser(option_list=option_list), args = t)
argv <- parse_args(OptionParser(option_list=option_list))#, args = t)
source("src/simulation_processing_tools.R")
source('../../../../src/plot_functions.R')
#source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
#source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/factorization_methods.R")
source("../../src/factorization_methods.R")
if(argv$ncores > 1)
{
  #Set up parallel run
  library(foreach)
  library(doParallel)
  n.cores=argv$ncores
  if(argv$ncores == -1)
  {
    n.cores <- parallel::detectCores() - 1
  }
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  #how many workers are available? (optional)
  foreach::getDoParWorkers()
  
  #Stop cluster when done
  #parallel::stopCluster(cl = my.cluster)
  
  
  
}
.libPaths("/home/aomdahl1/R/4.2.1")
print(.libPaths())
#assess via all the methods:
save.list <- list()
message("Starting file read in")
sim.path <- argv$input_dir
files.betas <- list.files(sim.path,pattern = "*.BETA.csv")
save.list$files.betas <- files.betas
sim.betas <- readInGWAS(files.betas, sim.path) 
rank <- ncol(sim.betas[[1]])
message("Read in a total of ", length(sim.betas), " files containing effect size estimates.")
files.se <- list.files(sim.path,pattern = "*\\.SE.csv")
sim.se <- readInGWAS(files.se, sim.path)

#Make sure the nubmers line up
stopifnot(length(sim.se) == length(sim.betas))
stopifnot(all(sapply(files.se, 
                     function(x) gsub(x, pattern = ".SE.csv", replacement = "")) == sapply(files.betas, 
                                                                                          function(x) gsub(x, pattern = ".BETA.csv", replacement = ""))))


sim.zscores <- lapply(1:length(sim.se), function(i) as.matrix(sim.betas[[i]])/as.matrix(sim.se[[i]]))

if(argv$run_only == "ALL" | grepl(pattern="PCA$", x=argv$run_only))
{
  save.list <- c(save.list, stdPCAEvals(sim.zscores))
}


if(argv$with_covar)
{
  message("Reading in covar and covar SE data")
  sim.cov <- readInGWAS(list.files(sim.path,pattern = "*\\.COV.csv"), sim.path) 
  sim.cov_se <- readInGWAS(list.files(sim.path,pattern = "*\\.COV_SE.csv"), sim.path) 
  stopifnot(length(sim.cov) == length(sim.cov_se))
  save.list$gleaner.runs.covar <- list()
}


if(argv$run_only == "ALL" | grepl(pattern="SVD_whiten", x=argv$run_only))
{
  sim.zscores.adj <- list()
  for(i in 1:length(sim.cov))
  {
    blocks <- create_blocks(sim.cov[[i]],cor_thr=0.2)
    covar <- blockifyCovarianceMatrix(blocks, as.matrix(sim.cov[[i]]))
    covar <- strimmerCovShrinkage(list(), covar, as.matrix(sim.cov_se[[i]]), sd.scaling=1)
    sim.zscores.adj[[i]] <- adjustMatrixAsNeeded(as.matrix(sim.betas[[i]])/as.matrix(sim.se[[i]]),covar)

  }
  save.list <- c(save.list, stdPCAEvals(sim.zscores.adj, prefix="whitened_"))
}

#Repeat all the above, except with 





message("Proceeding now with gleaner and flash")
save.list$flashr.runs.z <- list()
save.list$flashr.runs.se <- list()
save.list$gleaner.runs.std <- list()
if(argv$ncores == -1)
{
  if(argv$run_only == "ALL" | grepl(pattern="GLEANER$", x=argv$run_only) | grepl(pattern="FLASH$", x=argv$run_only))
  {
    for(i in argv$start:length(sim.betas))
    {
      message("On simulation number ",i )
      message(files.betas[[i]])
      if(argv$run_only == "ALL" | grepl(pattern="GLEANER$", x=argv$run_only))
      {
        suppressMessages(save.list$gleaner.runs.std[[i]] <- runSingle("GLEANER_glmnet_noCovar", as.matrix(sim.betas[[i]]), K="GRID", se_m=as.matrix(sim.se[[i]]), covar = NULL, bic.var = "sklearn_eBIC",
                                                                      init.mat = "V", is.sim = TRUE, 
                                                                      save.path = "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/sandbox/gleaner")) #not accounting for BIC typoe #Uhoh- we should be outputting something, otherwise the signal cutting is just too strong....
        #Works with Zou, but not with sklearn.... arg. Why why why why why why
        #Zou is a definite no-no here, I mean 8 factors is just too many.
        #Dev- zeeroes out
      }
      if(argv$run_only == "ALL" | grepl(pattern="FLASH$", x=argv$run_only))
      {
        suppressMessages(save.list$flashr.runs.z[[i]] <- flashr::flash(data=as.matrix(sim.betas[[i]])/as.matrix(sim.se[[i]]))) #uhoh!)
        ###Update 9/25- flash version with SE
        ###
        message("Std flash")
        suppressMessages(save.list$flashr.runs.se[[i]] <- flashier::flash(as.matrix(sim.betas[[i]]),
                                                                          S=as.matrix(sim.se[[i]]), backfit = TRUE,var_type = 0)) #2 empirically does poorer in simulations, so using 0#0 is constant- probably not best here...
       
        #message("flash + kronek")
        message("Omitting kronek for now")
        if(FALSE)
        {
          suppressMessages(save.list$flashr.runs.kr[[i]] <- tryCatch(
            {
              flashier::flash(as.matrix(sim.betas[[i]]),
                              S=as.matrix(sim.se[[i]]), backfit = TRUE,var_type = c(1,2))
            },
            error = function(cond) {
              message("Kronecker calc failed")
              NA
            },
            finally = {
              message("Kronecker calc done")
            }
          ))
        }
          
      
      }
      
      if(argv$with_covar)
      {
        if(argv$run_only == "ALL" | grepl(pattern="GLEANER$", x=argv$run_only))
        {
          suppressMessages(save.list$gleaner.runs.covar[[i]] <- runSingle("GLEANER_glmnet", as.matrix(sim.betas[[i]]), K="GRID", se_m=as.matrix(sim.se[[i]]), 
                                                                          covar = as.matrix(sim.cov[[i]]), covar_se=as.matrix(sim.cov_se[[i]]), bic.var = "dev",
                                                                          init.mat = "V", is.sim = TRUE, enforce_blocks=FALSE, shrinkWL="strimmer",conv_objective =argv$convergence_criteria,
                                                                          save.path = "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/sandbox/gleaner"))#not accounting for BIC typoe #Uhoh- we should be outputting something, otherwise the signal cutting is just too strong....
        }
      }
      if(i%%25 == 0)
      {
        #perform an intermediate save
        message("Intermediate update.")
        save(save.list,
             file=argv$output)
        
      }
    }
  }

  
}else
{
  message("running in parallel.")
  foreach(i = argv$start:length(sim.betas)) %dopar% 
  {
    print("On simulation number ",i )
    message(files.betas[[i]])
    
    if(argv$run_only == "ALL" | grepl(pattern="GLEANER$", x=argv$run_only))
    {
      suppressMessages(save.list$gleaner.runs.std[[i]] <- runSingle("GLEANER_glmnet_noCovar", as.matrix(sim.betas[[i]]), K="GRID", se_m=as.matrix(sim.se[[i]]), covar = NULL, bic.var = "dev",
                                                        init.mat = "V", is.sim = TRUE, 
                                                        save.path = "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/sandbox/gleaner")) #not accounting for BIC typoe #Uhoh- we should be outputting something, otherwise the signal cutting is just too strong....
    #Works with Zou, but not with sklearn.... arg. Why why why why why why
    #Zou is a definite no-no here, I mean 8 factors is just too many.
    #Dev- zeeroes out
      if(argv$with_covar)
      {
        suppressMessages(save.list$gleaner.runs.covar[[i]] <- runSingle("GLEANER_glmnet", as.matrix(sim.betas[[i]]), K="GRID", se_m=as.matrix(sim.se[[i]]), 
                                                                        covar = as.matrix(sim.cov[[i]]), covar_se=as.matrix(sim.cov_se[[i]]), bic.var = "dev",
                                                                        init.mat = "V", is.sim = TRUE, enforce_blocks=FALSE, shrinkWL="strimmer",conv_objective =argv$convergence_criteria,
                                                                        save.path = "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/sandbox/gleaner"))#not accounting for BIC typoe #Uhoh- we should be outputting something, otherwise the signal cutting is just too strong....
        
      }
      
    }
    if(argv$run_only == "ALL" | grepl(pattern="FLASH$", x=argv$run_only))
    {
      suppressMessages(save.list$flashr.runs[[i]] <- flashr::flash(data=as.matrix(sim.betas[[i]])/as.matrix(sim.se[[i]]))) #uhoh!)
    }

  }
  parallel::stopCluster(cl = my.cluster)
}

#Save the output data:
save(save.list,
     file=argv$output)
#save(flashr.runs,gleaner.runs, file="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/sandbox/flash_v_gleaner_N-var.RData")



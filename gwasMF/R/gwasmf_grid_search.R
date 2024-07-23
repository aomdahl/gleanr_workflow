#Source everything you need:

SourcePackages <- function()
{
  dir ="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/"
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
  source(paste0(dir, 'evaluateRunsMany.R'))
  source(paste0(dir, 'pve.R'))
}

quickSort <- function(tab, col = 1)
{
  tab[order(tab[,..col], decreasing = TRUE),]
}

updateStatement  <- function(l,a,l_og, a_og, run,time)
{

  #logr::log_print(paste0("Run complete. (Time: ", time, " min)"))
  #logr::log_print(paste0("Sparsity params: F ", l, ", L ", a))
  #logr::log_print(paste0("Sparsity scale: F ", l_og, ", L ", a_og))
  #logr::log_print(paste0("Current F sparsity: ", run$V_sparsity), console = FALSE)
  #logr::log_print(paste0("Current L sparsity: ", run$U_sparsity), console = FALSE)
  #logr::log_print(paste0("Number of active factors: ", ncol(run$V)))

  message(paste0("Run complete. (Time: ", time, " min)"))
  message(paste0("Sparsity params: F ", l, ", L ", a))
  message(paste0("Sparsity scale: F ", l_og, ", L ", a_og))
  message(paste0("Current F sparsity: ", run$V_sparsity), console = FALSE)
  message(paste0("Current L sparsity: ", run$U_sparsity), console = FALSE)
  message(paste0("Number of active factors: ", ncol(run$V)))
}


initializeGwasMF <- function(X,W,C,snp.ids, trait.names, K=0, init.mat = "V", covar_shrinkage=-1,enforce_blocks=TRUE,covar_se=NULL, ...)
{
  #SourcePackages()
  #A few quick sanity checks:
  if(all(apply(X, 2, var) == 0) | all(apply(X, 1, var) == 0))
  {
    message("Matrix of effect sizes contains no variance. Likely not a valid matrix- program will terminate")
    quit()
  }
  if(all(apply(W, 2, var) == 0) | all(apply(W, 1, var) == 0))
  {
    warning("Matrix of standard errors contains no variance. May not be valid")
  }

  message("This is an archaic initialization; recommend doing away with this...")
  args <- defaultSettings(K=K,init.mat = init.mat,...)
  args$pve_init <- FALSE
  option <- readInSettings(args)
  option$swap <- FALSE
  option$alpha1 <- 1e-10
  option$lambda1 <- 1e-10
  output <- args$output
  #log.path <- paste0(args$output, "gwasMF_log.", Sys.Date(), ".txt")
  #lf <- logr::log_open(log.path, show_notes = FALSE)
  #options("logr.compact" = TRUE)
  #options("logr.notes" = FALSE)
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

  #Adjusted blocking approach- block first, then shrinkge
  #this isn't needed in simulations, so I think we droop this
  if(enforce_blocks)
  {
    blocks <- create_blocks(C,cor_thr=0.2)
    covar <- blockifyCovarianceMatrix(blocks, C)
  }else
  {
    message("Not enforcing a block structure")
    covar <- C
  }

  #if specified, srhink
  #WHAT IS THE ORDER
  if(covar_shrinkage == "strimmer" | covar_shrinkage == "STRIMMER")
  {
    covar <- strimmerCovShrinkage(args, covar, covar_se, sd.scaling=1)
  }
  else if(covar_shrinkage > -1)
  {
    message("Performing WL shrinkage, as desired")

  covar <- linearShrinkLWSimple(covar, args$WLgamma)

  }
  else{
    message("Not applying shrinkage on covar")
  }
 W_c_dat <- buildWhiteningMatrix(covar, ncol(X),blockify = -1)

 # W_c_dat <- buildWhiteningMatrix(C, blockify = 0.2)
  W_c = W_c_dat$W_c
  #(list("W_c" = solve((chol(covar))),"C_block"=covar))
  list("args" = args, "options" = option,
       "hp" = hp, "all_ids" = snp.ids, "names" = trait.names, "initk" = initk, "W_c" = W_c)
}


#' Perform a single itertion of factorization with given options and arguments
#' I'm pretty sure this is redundant with other code, but here is again.
#'
#' @param X input data matrix
#' @param W uncertainy matrix
#' @param W_c decorrelation matrix
#' @param a alpha to run at
#' @param l lambda to run at
#' @param option option settings
#' @param niter how many iterations to allow.
#'
#' @return the full factorization return data
#' @export
singleRunFromScratch <- function(X,W,W_c,a, l, option, niter = 100)
{
  option[['alpha1']] <- as.numeric(a)
  option[['lambda1']] <- as.numeric(l)
  option$iter <- niter
  Update_FL(as.matrix(X), as.matrix(W),W_c, option)
}

#getGridMatrices(args, option, alphas, lambdas, param.iters)
#' Performs the grid expansion for alphas and lambdas introduced here..
#'
#' @param X data matrix
#' @param W uncertainty weightings
#' @param W_c Decorrelation matrix
#' @param args Argument list
#' @param option option list
#' @param alphas list of alphas to test
#' @param lambdas list of lambdas to test
#' @param param.iters Number of iterations to repeat each run at to get stability.
#'
#' @return list containing all the run data.
#' @export
getGridMatrices <- function(X,W,W_c,args, option, alphas, lambdas, param.iters)
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
        if(is.na(a) | is.na(l))
        {
          message("Proposed a, l are NA, skipping.")
          next;
        }
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
            X <- as.matrix(X)$X[samp,]; W <- as.matrix(W)$W[samp,]; #all_ids <- input.dat$ids[samp]
            #This works, but need to adjust the parameters. COuld implement okay, but not doing now.f
          }
        run_stats[[opti]] <- Update_FL(as.matrix(X), as.matrix(W), as.matrix(W_c),option)
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
#res <- gwasMFGrid(X,W, paste0("rs", 1:nrow(X)), paste0("T",1:ncol(X)), K=0, Kmax = K)
gwasMFGrid <- function(X,W, snp_ids, trait_names, C=NULL, K=0, Kmax = -1, subsample = FALSE, nrep.first = 5, nrep.second = 10)
{
  print("Setting for k")
  print(K)
  if(is.null(C))
  {
    C <- diag(ncol(X))
  }
  d <- initializeGwasMF(X,W, C,snp_ids, trait_names, K=0)
  option <- d$options; args <- d$args; hp <- d$hp; all_ids <- d$all_ids; names <- d$names; W_c <- d$W_c
  option$svd_init <- FALSE
  args$svd_init <- FALSE
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
  init.params <- GetCoarseSparsityParams(X,W,W_c,option, burn.in.reps = 2, logs = TRUE, n.points = 4) #TODO- find a faster way to do this.
  alphas <- init.params$alphas; lambdas <- init.params$lambdas
  message("Alphas")
  print(alphas)
  message("Lambdas")
  print(lambdas)

  grid.solutions <- getGridMatrices(X, W,W_c, args, option, alphas, lambdas, nrep.first)

  check_stats = max(sapply(grid.solutions, length))
  if(check_stats == 1)
  {
    #logr::log_print("No runs completed. Terminating")
    message("No runs completed. Terminating")
    quit()
  }

  ret <- GetNextParams(alphas, lambdas,init.k, grid.solutions)

  while(all(is.na(ret)))
  {
	  message("ERROR Error Got an NA case- this means the space we are searching is probably not valid?")
	  mesage("Trying step again, with less stringent parameters:")
	  alphas <- alphas * 0.1
	  lambdas <- lambdas * 0.1
	  grid.solutions <- getGridMatrices(X, W,W_c, args, option, alphas, lambdas, nrep.first)
	  ret <- GetNextParams(alphas, lambdas,init.k, grid.solutions)

  }
  option$conv0 <- 0.001
  grid.solutions.round.final <- getGridMatrices(X,W,W_c,args, option, ret$alphas, ret$lambdas, nrep.second)
  #Above is good.
  final.run <- GetGridPerformanceMetrics(ret$alphas, ret$lambdas, grid.solutions.round.final, init.k)

  best.params <- SelectTopRun(final.run, grid.solutions.round.final, maxK=Kmax)
  if(any(is.na(best.params)))
  {
    message("Unable to find final parameters.")
    return(NULL)
  }
  reg.run <- singleRunFromScratch(X,W,W_c, best.params$alpha, best.params$lambda, option, niter = 200)
  #function(X,W,W_c,U,V, minK, option, maxK = NULL,drop.min.change = TRUE)
  #d <- DropFactorsByObjective(X,W,W_c,reg.run$U,reg.run$V, minK=Kmax, option, maxK = Kmax)
  reg.run <- PruneNumberOfFactors(X,W,W_c,reg.run,Kmax, option$Kmin, option)

  #Question to answer- is it best to choose according to number of params and then parse it down, or how?
  #For now- pick the one that corresponds to your target number of factors
  #Above
  #return(list(d, reg.run))
  return(reg.run)
  }


GetGridPerformanceMetrics <- function(alphas, lambdas,grid.output, init.k)
{
  #if("dplyr" %in% (.packages())){
   # detach("package:dplyr", unload=TRUE)
  #  detach("package:plyr", unload=TRUE)
  #}
  #library(plyr)
 # library(dplyr)

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
  by.group.med <- final.out %>% dplyr::group_by(Group_nonempty_average) %>% summarise("Coph_med" = median(Coph))
  left_join(final.out, by.group.med, by = "Group_nonempty_average")
}


SelectTopRun <- function(perf.df, runs.by.name, maxK=-1)
{

  filtered.df <- perf.df %>% filter(Coph > 0.9, Coph_med > 0.9, Avg_pve > 0.5, Avg_pve < 1,
                     V_sparsity_raw > 0.05, U_sparsity_raw > 0.01 ) %>% arrange(Med_Average_R2)
  if(maxK != -1)
  {
    filtered.df <- filtered.df %>% filter(Group_nonempty_average <= maxK)

  }
  if(nrow(filtered.df) < 1)
  {
    filtered.df <- perf.df %>% filter(Coph > 0.9, Coph_med > 0.9)  %>% arrange(Med_Average_R2)
    if(nrow(filtered.df) < 1)
    {
      message("Unable to find stable meatrix")
      return(list("alpha"=NA, "lambda"=NA))
    }
  }
  #Return the top settings
  #Choosingt he one with the minium correlation
  return(list("alpha"=as.numeric(filtered.df[1,]$Alpha), "lambda"=as.numeric(filtered.df[1,]$Lambda)))
}
GetNextParams <- function(alphas, lambdas,init.k, all.solutions)
{
  perf.df <- GetGridPerformanceMetrics(alphas, lambdas,all.solutions, init.k)
  #print(perf.df)
  #save(perf.df, file = "debuggingGrid.RData")
  stopifnot(nrow(perf.df) > 0)
  #begin with most stringent selection criteria
  tops <- perf.df %>% tidyr::drop_na() %>% filter(Coph > 0.9) %>% filter(V_sparsity_global > 0.05, V_sparsity_raw > 0.05) %>%
             arrange(Med_Average_R2) %>% filter(!is.na(Coph_med), Coph_med > 0.9, Avg_pve > 0.5, U_sparsity_raw > 0.05)

  #how do we narrow this down to a small number of points to search?
  if(nrow(tops) >= 1)
  {
    best.a <- unique(as.numeric(gsub(tops$Alpha, pattern = "A", replacement = "")))
    best.l <- unique(as.numeric(gsub(tops$Lambda, pattern = "L", replacement = "")))
    if(length(best.a) > 2 & length(best.l) > 2)
    {
  	  message("Selection criteria not strict enough, parsing down further")
  	  tops <- (tops %>% filter( V_sparsity_raw > 0.1) %>% arrange(Med_Average_R2))[1:2,] #Just choosing the top, most sparse
    	best.a <- unique(as.numeric(gsub(tops$Alpha, pattern = "A", replacement = "")))
    	best.l <- unique(as.numeric(gsub(tops$Lambda, pattern = "L", replacement = "")))
    }
    new.a <- ProposeSparsityParamsFromGrid(best.a,alphas, 4)
    new.l <- ProposeSparsityParamsFromGrid(c(best.l),lambdas,4)
  }else #(nrow(tops) == 0)
  {
    message("Parameters too strict. Choosing the most stable solution with any sparsity and K > 1.")
    #Just choosing the most stable one with some sparsity
    tops <- data.frame(perf.df) %>% tidyr::drop_na() %>%
      filter(Group_nonempty_average > 1, V_sparsity_raw > 0) %>% arrange(-Coph,Med_Average_R2) %>% slice(1)
    if(nrow(tops) < 1) #still too strict
    {
	    print('things arent loooking good...')
	    print(perf.df)
	    tops <- data.frame(perf.df) %>% arrange(-Coph,Med_Average_R2) %>% slice(1)
    }
    best.a <- unique(as.numeric(gsub(tops$Alpha, pattern = "A", replacement = "")))
    best.l <- unique(as.numeric(gsub(tops$Lambda, pattern = "L", replacement = "")))
    new.a <- ProposeSparsityParamsFromGrid(best.a,alphas, 4)
    new.l <- ProposeSparsityParamsFromGrid(c(best.l),lambdas,4)
  }
  #dropNAs
  new.a <- new.a[!is.na(new.a)]
  new.l <- new.a[!is.na(new.l)]
  if(length(new.a) < 1 | length(new.l) < 1)
  {
    message("Error in parameter selection. Need to debug urgently.")
    message("Likely had NAs in int.")
    print(alphas)
    print(lambdas)
    return(NA)
  }
  return(list("alphas" = new.a, "lambdas" = new.l))

}


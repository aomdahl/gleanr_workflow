#Prep data

#scan the input
readInParamterSpace <- function(args)
{
  #Read in the hyperparameters to explore
  if(args$autofit)
  {
    log_print("Autofit setting of sparsity parameters detected. ")
    alphas_og <- c(NA)
    lambdas_og <- c(NA)
    if(args$alpha != "" | args$lambdas != "")
    {
      log_print("Incompatible sparsity settings provided- cannot autofit and use specified. Program will terminate.")
      quit()
    }
  }else
  {
    alphas_og <- as.numeric(scan(text = args$alphas, what = character(), sep = ',', quiet = TRUE))
    lambdas_og <- as.numeric(scan(text = args$lambdas, what = character(), sep = ',', quiet = TRUE))
  }

  return(list("a" = alphas_og, "l" = lambdas_og))
  
}

#Holdover from previous version. May not be used
cleanUp <- function(matin, type = "beta")
{
  print(dim(matin))
  lister <- as.vector(unlist(matin))
  if(type == "beta")
  {
    bad <- is.na(lister)
    lister[bad] <- 0
    bad <- is.infinite(lister)
    lister[bad] <- 0
  }
  if(type == "se")
  {
    bad <- is.infinite(lister)
    lister[bad] <- 1000
    bad <- is.na(lister)
    lister[bad] <- 1
  }
  return(matrix(lister, nrow = nrow(matin)))
}

#Zero out NA and extreme chi2
#Param X, W
#Return list with cleanX, cleanW
matrixGWASQC <- function(X, W)
{
  #Now that we have the SE and the B, go ahead and filter out bad snps....
  #Clean up NAs, etc.
  #Sanity check for number ofo NAs. X*W gives Z scores
  #CHECK  drops by SNP first
  tnames = colnames(X)
  drop.nas.rows <- apply(X*W, 1, function(x) sum(is.na(x)))
  drops = c()
  if(any(drop.nas.rows == ncol(X)))
  {
    X <- X[-which(drop.nas.rows == ncol(X)),]
    W <- W[-which(drop.nas.rows == ncol(X)),]
    drops = which(drop.nas.rows == ncol(X))
    log_print(paste0("Removed ",sum(drop.nas.rows == ncol(X)), " variants where all entries were NA..."))
  }
  drop.nas <- apply(X*W,2, function(x) is.na(x))
  drop.chi2 <- apply(X*W,2, function(x) x^2 > 80) 
  drop.chi2[drop.nas] <- FALSE
  all.drops <- drop.chi2 | drop.nas
  #do we just zero those ones out or drop them all together?
  #if more than 15% of SNPs for a particular trait are in this category, drop the trait instead
  too.many.drops <- unlist(lapply(1:ncol(drop.chi2), function(i) sum(drop.chi2[,i]) + sum(drop.nas[,i])))
  if(any(too.many.drops > 0.20 * nrow(X))) #greater than 20%
  {
    drop.counts <- which(too.many.drops > 0.20 * nrow(X))
    for(w in drop.counts)
    {
      
      log_print(paste0(round(too.many.drops[w]/nrow(X) * 100, digits = 2), "% of SNPs in trait ",tnames[w], " are either missing or invalid. Consider removing this trait."))
    }
  }
  log_print("Cells with invalid entries (NA, or excessive Chi^2 stat) will be given a weight of 0.")
  log_print(paste0(sum(all.drops), " out of ", (nrow(all.drops) * ncol(all.drops)), " total cells are invalid and will be 0'd."))
  log_print(paste0("   Zeroing out ", sum(drop.chi2), " entries with extreme chi^2 stats > 80"))
  log_print(paste0("   Zeroing out ", sum(drop.nas), " entries with NA"))
  W[all.drops] <- 0
  X[drop.nas] <- 0 #technically we don't need to drop them here, if they are being given a weight of 0. But if they are NAs we need to give them a value so they don't killus.
  if(length(drops) == 0)
  {
	  drops = 0
  }
  return(list("clean_X" = X, "clean_W" = W, "dropped_rows"=drops))
}



#Takes care of all the necessary data I/O
#Procedures include:
# Read in files
# Specify the weighting scheme
# "Drop" summary statistics with X^2 stat > 80 (set its weight to 0)
# Scale by sqrt(N), if specified
# Scale by sqrt(LDSC_int) if specified
#@param args- argument object in R
#@return a list containing the values, their corresponding weights, and the SNPS in the matching order
readInData <- function(args)
{
  #Load the effect size data
  effects <- fread(args$gwas_effects) %>% quickSort(.)
  all_ids <- unlist(effects[,1])
  
  #Get the trait names out
  if(args$trait_names == "")
  {
    message("No trait names provided. Using the identifiers in the tabular effect data instead.")
    names <- unlist(names(effects)[-1])
  } else{
    names <- scan(args$trait_names, what = character(), quiet = TRUE)
  }
  
  effects <- as.matrix(effects[,-1])
  #Do we IRNT normalize the input data
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
    W_se <- fread(args$uncertainty) %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.)
    stopifnot(all(all_ids == W_se[,1]))
    W_se <- W_se %>% select(-1)
    W <- 1/ W_se
    X <- effects
    
  } else if(args$weighting_scheme == "B_MAF")
  {
    message("Scaling by 1/var(MAF)")
    W_maf <- fread(args$uncertainty) %>%
      filter(ids %in% all_ids) %>% arrange(ids) %>% select(-ids) 
    W <- 1/matrix(apply(W_maf, 2, function(x) 2*x*(1-x)), nrow = nrow(W_maf), ncol = ncol(W_maf))
    X <- effects 
  } else
  {
    message("No form selected. Please try again.")
    quit()
  }
  
  #remove NAs, extreme values.
  r <- matrixGWASQC(X,W)
  X <- r$clean_X;   W <- r$clean_W; all_ids <- all_ids[-(r$dropped_rows)]
  #browser()
  #At this point, we are assuming no ambiguous SNPs, but maybe we shouldn't....
  #NO AMBIGS
  #NO MAF < 0.01
  
  if(args$scale_n != "")
  {
    N <- fread(args$scale_n) %>% filter(!row_number() %in% r$drops) %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.)
    if(!all(all_ids == N[,1]))
    {
      #This would occur if we are missing variants.
      message("Counts not provided for all variants. Using the average where variants missing")
      vars <- all_ids[!(all_ids %in% unlist(N[,1]))]
      m <- colMeans(N[,-1])
      ndiff <- nrow(X) - nrow(N)
      pre <- lapply(1:ndiff, function(x) unlist(m))
      
      first <- do.call("rbind", pre)
      new_rows <- cbind("SNP" = vars,first)
      N <- rbind(N, new_rows) %>% quickSort(.)
      stopifnot(all(all_ids == N[,1]))
    }
    N <- as.matrix(N %>% select(-1) %>% apply(., 2, as.numeric))
    W <- W * (1/sqrt(N))
  }
  
  if(args$genomic_correction != "")
  {
    log_print("Including genomic correction in factorization...")
    GC <-  fread(args$genomic_correction) %>% filter(!row_number() %in% r$drops) %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.)
    if(!all(all_ids == GC[,1])) {
      message("genomic correction values not in correct order, please advise...")
      GC <- as.matrix(GC %>% select(-1) %>% apply(., 2, as.numeric))
      W <- W * (1/sqrt(GC))
    }
  }
  

  return(list("X" = X, "W" = W, "ids" = all_ids, "trait_names" = names))
  
}


readInSettings <- function(args)
{
  option <- list()
  option[['K']] <- args$nfactors
  option[['iter']] <- args$niter
  option[['convF']] <- 0
 option[["nsplits"]] <- as.numeric(args$cores)
  option[['convO']] <- args$converged_obj_change
  option[['ones']] <- FALSE
  option[['disp']] <- FALSE
  #F matrix initialization
  option[['f_init']] <- args$init_F
  option[['epsilon']] <- as.numeric(args$epsilon)
  option[['l_init']] <- args$init_L
  option[["preinitialize"]] <- FALSE
  option[['reweighted']] <- FALSE
  option[["glmnet"]] <- FALSE
  option[["parallel"]] <- FALSE
  option[["fastReg"]] <- FALSE
  option[["ridge_L"]] <- FALSE
  
  if(args$regression_method == "penalized" | args$regression_method == "glmnet")
  {
	 option[["regression_method"]] = args$regression_method #push this through all initializations.
  }else{
	  message("This method isn't recognized. Try penalized or glmnet")
	  quit()
  }
  option[["posF"]] <- args$posF
  option[["autofit"]] <- args$autofit
  option$intercept_ubiq <- FALSE
  option$traitSpecificVar <- FALSE
  option$V <- args$verbosity
  option$calibrate_sparsity <- args$scaled_sparsity
  if(args$cores > 1)
  {
    log_print(paste("Running in parallel on", args$cores, "cores"))
    option[["parallel"]] <- TRUE
  }
  option[["ncores"]] <- args$cores
  option[["fixed_ubiq"]] <- args$fixed_first
  return(option)
}

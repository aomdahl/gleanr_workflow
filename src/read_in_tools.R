#Prep data

#scan the input
readInParamterSpace <- function(args)
{
  #Read in the hyperparameters to explore
  alphas_og <- as.numeric(scan(text = args$alphas, what = character(), sep = ',', quiet = TRUE))
  lambdas_og <- as.numeric(scan(text = args$lambdas, what = character(), sep = ',', quiet = TRUE))
  return(list("a" = alphas_og, "l" = lambdas_og))
  
}


readInData <- function(args)
{
  #Load the effect size data
  effects <- fread(args$gwas_effects) %>% drop_na() %>% quickSort(.)
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
    W_se <- fread(args$uncertainty) %>% drop_na() %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.)
    stopifnot(all(all_ids == W_se[,1]))
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
  
  #At this point, we are assuming no ambiguous SNPs, but maybe we shouldn't....
  #NO AMBIGS
  #NO MAF < 0.01
  #Now that we have the SE and the B, go ahead and filter out bad snps....
  #Clean up NAs, etc.
  #Sanity check for number ofo NAs
  drop.chi2 <- apply(X*W,2, function(x) which(x^2 > 80))
  drop.nas <- apply(X*W,2, function(x) which(is.na(x)))
  #do we just zero those ones out or drop them all together?
  #if more than 10% of SNPs for a particular trait are in this category, drop the trait instead
  which(lapply(drop, function(x) length(drop))  > 0.1 * nrow(X))
  
  drop <- (apply(b[,-1]/se[,-1], 2, function(x) which(x^2 > 80)))
  drop.traits <- which(lapply(drop, function(x) length(drop))  > 0.1 * nrow(z))
  drop.vars <- unique(unlist(drop))
  if(length(drop.vars) > ( 0.1 * nrow(z)))
  {
    message("Over 10% of SNPs have highly unusual Z-scores. WE recommend reveiwing your data.")
    message("Rather than dropping them all, we will simply zero-out those particular entries.")
  }
  
  
  log_print(paste0("Removing ", length(drop), " variants with extreme chi^2 stats > 80"))
  X <- X[-drop,]
  W <- W[,-drop]
  all_ids <- all_ids[-drop]
  
  if(args$scale_n != "")
  {
    N <- fread(args$scale_n) %>% drop_na() %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.)
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
    GC <-  fread(args$genomic_correction) %>% drop_na() %>% filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.)
    if(!all(all_ids == GC[,1])) {
      message("genomic correction values not in correct order, please advise...")
      GC <- as.matrix(GC %>% select(-1) %>% apply(., 2, as.numeric))
      W <- W * (1/sqrt(GC))
    }
  }
  
  #max_sparsity <- Update_FL(as.matrix(X), as.matrix(W), option)
  X <- cleanUp(X)
  W <- cleanUp(W, type = "se")
  
  return(list("X" = X, "W" = W, "ids" = all_ids))
  
}

#Projection helper functions
#These come from Yuan He's code, and have been modified for my use here....

#Weight Some matrix by W. Flexible to take either F or X
#X- the matrix to weight. Can be the same dimensions of W, or just share one dimension
#W - the weights.
#z_score: if you want to force it to just do z scores
#Decorrelate- a covariance structure to adjust for by whitening. If NULL, no adjustment will occur
#precalc.U- the matrix to whiten by- this is an option if you've already calculated U, to avoid doing so again.
#weightEffects(factors, W, decorrelate,precalc.U = whitening.matrix )
#weightEffects(projection.target, W,decorrelate, precalc.U = whitening.matrix) 
weightEffects <- function(X,W, decorrelate, z_score = FALSE, precalc.U = NULL)
{
  source("/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/factorization_methods.R")
  #We presume if W is Ses, that W <- 1/SE
  #X cam 
  #Return item as a list
  if(is.null(W))
  {
    message("WARNING: null W")
    return(lapply(1:nrow(X), function(i) X[i,]))
           
  } else if(all(dim(X) == dim(W)) | z_score) {
    #return(lapply(1:nrow(X), function(i) X[i,]*W[i,]))
    #this takes too long, doing it another way
    #return(split(X*W,seq(nrow(X))))
    return(adjustMatrixAsNeeded(X*W, decorrelate, whitener = precalc.U)) #this is a matrix. if decorrelate is NULL, just returns the input.
  } else { 
    N <- nrow(W)
    #Code nippet from YUAn
    weighted_data = list()
      for(n in seq(1,N)){
        weighted_data[[n]] = list()
        w = as.numeric(W[n, ])
        weighted_data[[n]] = adjustMatrixAsNeeded(diag(w) %*% X, decorrelate, whitener = precalc.U)
      }
    return(weighted_data)
  }
}

#This function from YUAN, with tweaks
#Read in a factor matrix, and return a vector of trait names, and the composition per factor.
readInFactors <- function(input_file, scale, trait_names=NA,thresh = 0.01){
  ### read in the factor matrix
  ## if there are multiple runs, use the one whose factors being the most independent
  #OMitted that paart
  factors <- fread(input_file)
  traits <- factors[,1]
  
  #Make sure the names are right
  if(any(is.na(trait_names)))
  {
    message("Assuming traits in correct order")
  }else
  {
    #ensure the factor names and rownames are in the right order....
    if(!all(trait_names == factors[,1]))
    {
      factors$n <- factor(unlist(factors[,1]), levels = trait_names)
      factors <- factors %>% arrange(n) %>% select(-n)
      traits <- factors[,1]
    }
    #Omitted- names not included on the factors file anymore. Should do though.
    #after fixing, check again
    name_check <- any(factors[,1] != trait_names)
    if(name_check) {
      print("Error in names; please check input data")
      print(factors[,1])
      print(names(combined))
      print("Error in names")
      return(NA)
      #quit(save = "no")
    }
  }
  FactorM_optimal <- factors[,-1]
  
  ## scale the factors to have the same magnitude (same max abs value)
  if(scale)
  {
    scale_max = median(apply(abs(FactorM_optimal), 2, max))
    FactorM_optimal = apply(abs(FactorM_optimal), 2, function(x) scale_max / max(x) * x) * sign(FactorM_optimal)
  }

  K = ncol(FactorM_optimal)
  
    ## obtain the representative tissues/features for each factor
  #tissues = rownames(FactorM_optimal)
  #now traits
  compositions = apply(FactorM_optimal, 2, function(x) unlist(traits[which(x>thresh)]))
  
  return(list("F"=FactorM_optimal, "composition"=compositions, "trait"=unlist(traits)))
}
#name_order <- factor.dat$trait
readInCovariance <- function(p, name_order)
{
  if(p == "") {return(NULL)}
  else
  {
    message("We make the strong assumption that the matrix is in the correct order.")
    
    w.in <- fread(p) 
    pick.names <- which(colnames(w.in) %in% name_order)
    as.matrix(w.in[pick.names, ..pick.names])
  }
  

}
#Now actually do the projection
#This is analagous to what Yuan does in her method, pusing the learned factor matrix against the full set of SNPs
#TODO: allow for a whitening matrix here to adjust for overlapping effects. I think this might be the solution, if its needed.

#projected_data <- lmProjection(as.matrix(factors), combined, NULL, factor.dat$trait, factor.dat$composition, decorrelate = whiten.dat)
lmProjection <- function(factors, projection.target, W, trait.list, f.composition, decorrelate = NULL, make.blocks = FALSE)
{
  source("/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/factorization_methods.R")
  #toy W:
  if(is.null(W))
  {
    W <- matrix(rep(1, nrow(projection.target)* ncol(projection.target)), 
                nrow = nrow(projection.target), ncol = ncol((projection.target)))
  }
  #Get the matrix for whitening, if it exists
  whitening.matrix <- buildWhiteningMatrix(decorrelate, blockify = make.blocks)
  #^We do this here so we don't have to repeat each time on a step that is already SLOW.
  message("Weighting F")
  weighted.F <- weightEffects(factors, W, decorrelate,precalc.U = whitening.matrix ) #This returns a list
  message("Weighting X")
  weighted.X <- weightEffects(projection.target, W,decorrelate, precalc.U = whitening.matrix) #this returns a matrix
  message("Prepping object")
  #process.obj <- lapply(1:nrow(weighted.X), function(i) list("X" = weighted.X[i,], "F" = weighted.F[[i]]))
  return(run_regression(weighted.X, weighted.F, trait.list, f.composition))
}

#Taken directly from YUAN
#  process.obj <- lapply(1:nrow(weighted.X), function(i) list("X" = weighted.X[i,], "F" = weighted.F[[i]]))
#return(run_regression(process.obj, trait.list, f.composition))
# factor.dat$trait, factor.dat$composition
run_regression <- function(weighted.X, weighted.F, tissues, compositions){
  options(warn=-1)
  message("Performing regression...")
  library(parallel)
  #runid = paste0(Sys.time()) %>% gsub(., pattern = "[ :-]", replacement = "")
  if(nrow(weighted.X) > 100000) #the data is pretty big.
  {
    #doing it blockwise to save progress along the way, allow for stops.
    #Mostly a practical consideration.
    lim = 10000
	  Betas = matrix(NA, nrow = (lim), ncol = ncol(weighted.F[[1]]))
	  pValues = matrix(NA, nrow = (lim), ncol = ncol(weighted.F[[1]]))
	  se = matrix(NA, nrow = lim, ncol = ncol(weighted.F[[1]]))
	  message("Data is large. Going to perform blockwise")
	  start.point <- checkForExistingRuns()
	  n.snps <- nrow(weighted.X)
	  if( start.point < n.snps)
	  {
	    for(i in start.point:n.snps)
	    {  #res = compute_p(weighted_data[[i]], tissues, compositions)
	      res = compute_p( list("X" = weighted.X[i,], "F" = weighted.F[[i]]), tissues, compositions) #Check this.
	      j = ifelse(i > lim, i %% lim, i) %>% ifelse(. == 0, 2000, .)
	      Betas[j,] = res$B
	      pValues[j,] =res$P
	      se[j,] = res$SE
	      if(i %% lim == 0 | i == n.snps)
	      {
	        message("Reached iteration: ",i )
	        write.table(data.frame(round(Betas, digits = 5)), file = paste0("betas_", i, ".tmp.csv"),sep = ",", quote = FALSE, row.names = FALSE)
	        write.table(data.frame(round(se,digits = 5)), file = paste0("se_", i, ".tmp.csv"),sep = ",", quote = FALSE, row.names = FALSE)
	        write.table(data.frame(pValues), file = paste0("pvals_", i, ".tmp.csv"),sep = ",", quote = FALSE, row.names = FALSE)
	        Betas = matrix(NA, nrow = lim, ncol = ncol(weighted.F[[1]]))
	  	   se = matrix(NA, nrow = lim, ncol = ncol(weighted.F[[1]]))
	        pValues = matrix(NA, nrow = lim, ncol = ncol(weighted.F[[1]]))
	      }
	    }
	  }
    #Now, we've finished the runs.
	  f <- list.files(path = "./", pattern = "*.tmp.csv")
	  #we have some excess rows that are just NAs, remove those
	  Betas <- joinFiles(f, "betas_")[1:n.snps,]
	  write.table(data.frame(round(Betas, digits = 5)), file = paste0("betas_", i, ".complete.csv"),sep = ",", quote = FALSE, row.names = FALSE)
	  SE <- joinFiles(f, "se_")[1:n.snps,]
	  write.table(data.frame(SE), file = paste0("se_", i, ".complete.csv"),sep = ",", quote = FALSE, row.names = FALSE)
	  pValues <- joinFiles(f, "pvals_")[1:n.snps,]
	  write.table(data.frame(pValues), file = paste0("pvals_", i, ".complete.csv"),sep = ",", quote = FALSE, row.names = FALSE)
	  #remove all the files..
	  cleanUpFiles(f)
  }else
  {
	 #result = mclapply(weighted_data, function(x) compute_p(x,tissues, compositions),mc.preschedule = TRUE)
  	result = lapply(weighted_data, function(x) compute_p(x, tissues, compositions))
  	pValues = c()
  	Betas   = c()
  	SE = c()
  	for(r in result){
  	  Betas = rbind(Betas, r[[1]])
  	  pValues = rbind(pValues, r[[3]])
  	  SE = rbind(SE, r[[2]])
  	}
  }
  Betas = as.data.frame(Betas)
  pValues = as.data.frame(pValues)
  SE = as.data.frame(SE)
  #SE = NULL
  
  return(list("B"=Betas,"SE" = SE, "P" =pValues))
}

#Taken, with some tweaks
#compute_p( list("X" = weighted.X[i,], "F" = weighted.F[[i]]), tissues, compositions)
#I wonder if I don't want to quite follow this, because I don't want to induce so much sparsity.
#I want it to be nomrally distributed, for LDSC.
compute_p <- function(dataPoint, tissues, compositions){
  thresh = 1e-5
  #### run linear regression
  factor_exist = seq(1, dim(dataPoint[[2]])[2])
  invalid_idx = 0
  # remove factors where no tissue has non-NA (0) values in it, and assign 0 to tissues with NA value
  if(sum((dataPoint[[1]] == 0)) > 0){
    valid_idx = which((dataPoint[[1]] != 0))
    invalid_idx = which((dataPoint[[1]] == 0))
    data_exist = tissues[valid_idx]
    #This tells us which factors have at least 1/2 of their traits present (i.e. non-zero) in the current SNP data.
    factor_exist = which(sapply(compositions, function(x) length(intersect(data_exist, unlist(x)))>length(x)/2))

    dataPoint[[1]][invalid_idx] = 0
    #just in case the data direction got swapped
    if(length(dataPoint[[1]]) != length(tissues)){
      dataPoint[[1]] = t(dataPoint[[1]])
    }
    dataPoint[[2]][invalid_idx,] = 0 #0 out this trait, we want no contribution from it since it has a 0 effect here.
    dataPoint[[2]] = dataPoint[[2]][, factor_exist]#select out only traits that are present in the factor
    if(length(dataPoint[[2]]) == length(tissues)){
      dataPoint[[2]] = t(t(dataPoint[[2]]))
    }

    
    # remove the shared factor if there are other factors become co-linear with it (because of NA values)/other factors that are also fully loaded on every trait.
    #modified so its not strictly 0, but some non-zero elements allowed (to value of 1e-5)
    if(FALSE)
    {
      #decided to remove this. Don't want this top happen.
      if((1 %in% factor_exist) & (length(factor_exist) > 1)){
        temp_d2 = as.matrix(dataPoint[[2]][, 2:length(factor_exist)])
        #If the number of non-zero entries in factor 1   == the number of non-zero entries in any other factor
        if(sum(abs(dataPoint[[2]][, 1]) > thresh) == max(apply(temp_d2, 2, function(x) sum(abs(x) > thresh)))){
          #Drop that factor
          dataPoint[[2]] = dataPoint[[2]][, seq(2, length(factor_exist))]
          factor_exist = factor_exist[factor_exist!=1]
        }
      }
    }

    
    if(length(factor_exist) == 0){
      message("Empty factor case...")
      beta = rep(0, length(compositions))
      pv = rep(-1, length(compositions))
      se = rep(-1, length(compositions))
      #return(list(beta, pv))
      return(list("B"=as.numeric(beta), "SE"=as.numeric(se), "P"=as.numeric(pv)))
    }
  }
  
#dropped the shared sign requirement.
  
  
  lmfit = lm(unlist(dataPoint[[1]])~0+dataPoint[[2]], na.action = na.omit)
  beta = unlist(coef(lmfit))
  se = unlist(summary(lmfit)$coef[,2])
  pv = coef(summary(lmfit))[,4]
  
  if(length(pv) < length(beta)){
    pv[names(beta)[which(is.na(beta))]] = -1
    names(pv) = matrix(unlist(strsplit(names(pv), "]]")), nrow=length(pv), byrow = T)[,2]
    pv = pv[order(as.numeric(names(pv)))]
  }
  
  beta[is.na(beta)] = 0
  pv[is.na(pv)] = -1
  se[is.na(se)] = -1
  if(length(pv) < length(compositions)){
    beta_formatted = rep(0, length(compositions))
    pv_formatted = rep(-1, length(compositions))
    se_formatted = rep(-1, length(compositions))
    beta_formatted[factor_exist] = beta
    pv_formatted[factor_exist] = pv
    se_formatted[factor_exist] = se
    beta = beta_formatted; pv = pv_formatted; se = se_formatted
  }
  return(list("B"=as.numeric(beta), "SE"=as.numeric(se), "P"=as.numeric(pv)))
}


checkForExistingRuns <- function()
{
  f <- list.files(path = "./", pattern = "*.tmp.csv")
  if(length(f) == 0)
  {
    return(1)
  }
  else
  {
    return(max(sapply(f, function(x) as.numeric(str_match(x, "_([0-9]+).tmp")[2]))) + 1)
  }
}


joinFiles <- function(f, str)
{
  
  k <- which(sapply(f, function(x) grepl(x, pattern = str)))
  #Note- its important that we preserve order here....
  #Put the names in the right order now!
  nums <-  sapply(f[k] ,function(x) as.numeric(str_match(x, "_([0-9]+).tmp")[2]))
  new.order <- nums[order(nums)]
  #stopifnot(all(order(nums) == 1:length(nums)))
  stopifnot(all(order(sapply(names(new.order), function(x) as.numeric(str_match(x, "_([0-9]+).tmp")[2]))) == 1:length(nums)))
  do.call("rbind", lapply(names(new.order), function(x) fread(x)))

}

cleanUpFiles <- function(f)
{
  lapply(f, unlink)
}


#sanity check
dat <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/infertility/p0.001_FULL/flash_backfit_zscores/LM_projection/projected_hapmap3_loadings.txt")

#Helper functions for reading in and parseing code.
#TODO: develop some testing cases for each of these
quickSort <- function(tab, option, col = 1)
{
  if(!option$sort)
  {
    return(tab)
  }
  #tab <- as.matrix(tab)
  tab[order(tab[,col, with = FALSE], decreasing = TRUE),] #the package version.
  #tab[order(tab[,..col], decreasing = TRUE),] #the debug version
}


readInParamterSpace <- function(args)
{
  #Read in the hyperparameters to explore
  if(args$MAP_autofit > -1 | args$auto_grid_search)
  {
    updateLog("Sparsity parameters will be proposed by software.")
    alphas_og <- c(NA)
    lambdas_og <- c(NA)
    if(args$alpha != "" | args$lambdas != "")
    {
      updateLog("Incompatible sparsity settings provided- cannot autofit and use specified. Program will terminate.")
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
matrixGWASQC <- function(X, W, id.list, na_thres = 0.5, is.sim = FALSE)
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
    drops = which(drop.nas.rows == ncol(X))
    X <- X[-drops,]
    W <- W[-drops,]
    id.list <- id.list[-drops]

    updateLog(paste0("Removed ",sum(drop.nas.rows == ncol(X)), " variants where all entries were NA..."))
  }
  X <- as.matrix(X)
  W <- as.matrix(W)
  drop.nas <- apply(X*W,2, function(x) is.na(x))
  chi.thresh = 300
  drop.chi2 <- apply(X*W,2, function(x) x^2 > chi.thresh)
  drop.chi2[drop.nas] <- FALSE
  all.drops <- drop.chi2 | drop.nas
  #do we just zero those ones out or drop them all together?
  #if more than 15% of SNPs for a particular trait are in this category, drop the trait instead
  too.many.drops <- unlist(lapply(1:ncol(drop.chi2), function(i) sum(drop.chi2[,i]) + sum(drop.nas[,i])))
  drop.cols <- c()
  if(any(too.many.drops > 0.20 * nrow(X))) #greater than 20%
  {
    warning.counts <- which(too.many.drops > (0.20 * nrow(X)))
    drop.cols <- which(too.many.drops > (na_thres * nrow(X)))
    for(w in warning.counts)
    {
      na.proportion <- too.many.drops[w]/nrow(X)
      updateLog(paste0(round(too.many.drops[w]/nrow(X) * 100, digits = 2), "% of SNPs in trait ",tnames[w], " are either missing or invalid."))
    }
    updateLog(paste0("Traits with over ", na_thres*100, "% missing entries will be dropped automatically."))
  }

  #TODO: fix this report, it isn't quite right.
  updateLog("Cells with invalid entries (NA, or excessive Chi^2 stat) will be given a weight of 0.")
  removed.cols <- length(drop.cols) * nrow(X)
  updateLog(paste0(sum(all.drops), " out of ", (nrow(all.drops) * ncol(all.drops)), " total cells are invalid and will be 0'd."))
  updateLog(paste0("   Zeroing out ", sum(drop.chi2), " entries with extreme chi^2 stats > ", chi.thresh))
  updateLog(paste0("   Zeroing out ", sum(drop.nas), " entries with NA"))
  W[all.drops] <- 0
  X[drop.nas] <- 0 #technically we don't need to drop them here, if they are being given a weight of 0. But if they are NAs we need to give them a value so they don't killus.
  if(length(drops) == 0){	  drops = 0  }
  if(length(drop.cols) == 0){
    return(list("clean_X" = X, "clean_W" = W, "dropped_rows"=drops, "dropped_cols" =  0, "snp.list"=id.list))
    }
    return(list("clean_X" = X[, -drop.cols], "clean_W" = W[, -drop.cols], "dropped_rows"=drops, "dropped_cols" =  drop.cols, "snp.list"=id.list))

}



#tidy up the column-names. fread has some wierd behavior I don't like....

#Takes care of all the necessary data I/O
#Procedures include:
# Read in files
# Specify the weighting scheme
# "Drop" summary statistics with X^2 stat > 80 (set its weight to 0)
# Scale by sqrt(N), if specified
# Scale by sqrt(LDSC_int) if specified
#@param args- argument object in R
#@return a list containing the values, their corresponding weights, and the SNPS in the matching order
#readInBetas(args$gwas_effects, args)
readInBetas <- function(fpath, option)
{
  `%>%` <- magrittr::`%>%`
  data.table::fread(fpath, check.names = TRUE) %>% quickSort(.,option)
}

#Function to read in a covariance matrix.
#This code was copied from "projection_regression_helper.R", but belongs here
# readInCovariance(args$covar_matrix, names)
readInCovariance <- function(p, name_order)
{
  if(p == "" | is.null(p)) {return(NULL)}
  else
  {
    message("We make the strong assumption that the matrix rows and columns are in the same order. (check by looking at diags)")
    message("Enforcing the column order to be the same as input file name order")

    w.in <- as.matrix(data.table::fread(p, check.names = TRUE))
    row.names(w.in) <- colnames(w.in)
    stopifnot(all(diag(as.matrix(w.in)) == 1))
    if(!isSymmetric(w.in, tol = 1e-3))
    {
      #could just be due to numerical errors, round it
      if(!isSymmetric(round(w.in, digits = 2)))
      {
        message("WARNING: matrix may not be symmetric. Verify input covariance structure")
      }
    }
    #pick.names <- which(colnames(w.in) %in% name_order)
    #need to match the names
    as.matrix(w.in[name_order, name_order])
  }


}


readInLambdaGC <- function(fpath,X, names)
{
  `%>%` <- magrittr::`%>%`
  #note- this can come in as a matrix or as a simple list
  GC <-  data.table::fread(fpath, check.names = TRUE)
  if(nrow(GC) >= ncol(X)) #we have a single entry for each
  {
    message("Testing that the ordering is the same")
    if(nrow(GC) > ncol(X))
    {
      print(names)
      GC <- dplyr::filter(GC, phenotype %in% names)
      if(nrow(GC) != length(names))
      {
        missing <- names[which(!(names %in% GC$phenotype))]
        print("we are missing")
        print(missing)
      }
    }
    ecol = 1
    if(ncol(GC) == 2)  {   ecol = 2 } #the first one is names
    if(!all(unlist(GC[,1]) == names))
    {
      if(any(!(unlist(GC[,1]) %in% names)))
      {
        message("Mismatch in names, some don't align. Please check this")
        return(NA)
      }
      GC <- GC %>% dplyr::mutate("pheno_rank" = factor(phenotype, levels = names)) %>% arrange(pheno_rank) %>% select(-pheno_rank)
      stopifnot(GC$phenotype == names)
    }
    lambdagc <- as.matrix((do.call("rbind", lapply(1:nrow(X), function(x) unlist(GC[,ecol])))), nrow = nrow(X))
    GC <- lambdagc
  } else { #matrix version.
    GC <- GC %>% dplyr::filter(!row_number() %in% r$drops) %>% dplyr::filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.,args)
    if(!all(all_ids == GC[,1])) {
      message("genomic correction values not in correct order, please advise...")
      GC <- as.matrix(GC %>% select(-1) %>% apply(., 2, as.numeric))
      GC <- GC[,-r$dropped_cols]
    }

  }
  GC
}

SpecifyWeightingScheme <- function(effects, all_ids,all_phenos, args)
{
  effects <- as.matrix(effects[,-1]) %>% orderColumnsByName(., all_phenos,force.ref = args$trait_names)
  #Look at the weighting scheme options...
  if(args$weighting_scheme == "Z" || args$weighting_scheme == "B")
  {
    message("No scaling by standard error will take place. Input to uncertainty being ignored.")
    W <- matrix(1,nrow = nrow(effects), ncol = ncol(effects))
    X <- effects

  } else if(args$weighting_scheme == "B_SE")
  {
    message("Scaling by 1/SE.")
    W_se <- data.table::fread(args$uncertainty, check.names = TRUE) %>%
      dplyr::filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.,args)
    stopifnot(all(all_ids == W_se[,1]))
    W_se <- W_se[,-1] %>% orderColumnsByName(., all_phenos, force.ref = args$trait_names)
    message("Dims of W_se are now:", dim(W_se))
    W <- 1/ W_se
    X <- effects

  } else if(args$weighting_scheme == "B_MAF")
  {
    message("Scaling by 1/var(MAF)")
    W_maf <- data.table::fread(args$uncertainty, check.names = TRUE) %>%
      dplyr::filter(ids %in% all_ids) %>% arrange(ids) %>% select(-ids)
    W <- 1/matrix(apply(W_maf, 2, function(x) 2*x*(1-x)), nrow = nrow(W_maf), ncol = ncol(W_maf))
    X <- effects
  } else
  {
    message("No form selected. Please try again.")
    quit()
  }
  return(list("X" = X, "W" = W))
}


#' Function to order the columns by matching phenotype names
#'.
#' @param query.mat the matrix to check the order on
#' @param ref the reference order of phenotypes
#' @param force.ref the matrix doesn't yet have names, so just force the names in the order given in ref
#'
#' @return an ordered matrix
#' @export
orderColumnsByName <- function(query.mat, ref, force.ref = "")
{
  ret.mat <- as.matrix(query.mat) #matrix to return with right names
  if(force.ref == "") #if we learned the names on the fly from the effects file. Assumes the names ARE present in the file..
  {
    ret.mat <- ret.mat[,setColOrder(colnames(query.mat), ref)]
  }

  colnames(ret.mat) <- ref
  ret.mat
  #return(query.mat[,..o])
}

#' Title
#'
#' @param query
#' @param ref
#'
#' @return
#' @export
#'
#' @examples
setColOrder <- function(query, ref)
{
  #The reference might not include all the phenotypes, huh?
  ret.list <- 1:length(query)
  #best case- same entries
  if(length(query) == length(ref) & all(query == ref))
  {
    return(ret.list)
  }else if(all(ref %in% query)) #next best- its all there, just need to rearrange
  {
    reorder.q <- order(factor(query, levels = ref))
    return(reorder.q)
  }else #bad case
  {
    message("Phenotype names in query file don't match the reference phenotypes. Check this")
    print(query)
    print("")
    print(ref)
    quit()
  }
}
#' Load the relevant datasets for analysis based on an argument vector
#'
#' @param args Specifies the paths to files and program settings to run. Includes *gwas_effects, weighting_scheme, uncertainty, genomic_correction,* and *covar_matrix arguments*
#'
#' @return A list with entries:
#' "X": sorted GWAS effect sizes,
#' "W": sorted GWAS uncertainty weights (1/SE),
#' "ids": sorted ID order, "trait_names" = names of the traits,
#' "C": the read in covariance matrix and "W_c": the (blockified) matrix to be used for decorrelation
#' @export
readInData <- function(args)
{
  `%>%` <- magrittr::`%>%`
  #Load the effect size data
  rg <- NULL
  effects <- readInBetas(args$gwas_effects, args)

    all_ids <- unlist(effects[,1])
    #Get the trait names out
  if(args$trait_names == "")
  {
    message("No trait names provided. Using the identifiers in the tabular effect data instead.")
    names <- unlist(names(effects)[-1]) %>% make.names(.)
  } else{
    message("Using the provided trait names, and assuming all files have columns in the correct order.")
    message("It is your responsibility to verify this.")
    names <- scan(args$trait_names, what = character(), quiet = TRUE) %>% make.names(.)
  }
  if(length(names) > (ncol(effects) - 1))
  {
    message("Error- passed in effect size file and list of phenotype names are of different lengths.")
    quit()
  }
  weighted.dat <- SpecifyWeightingScheme(effects, all_ids,names, args)
  X <- weighted.dat$X; W <- as.matrix(weighted.dat$W);
  if(args$rg_ref != "")
  {
    message("Using an LDSC-rg based V initialization")
    rg <- readInCovariance(args$rg_ref, names)
  }
  if(args$drop_phenos != "")
  {
    message("Dropping data corresponding to names: ", args$drop_phenos)
    drops <- unlist(strsplit(args$drop_phenos, split = ",")[[1]])
    drop.indices <- which(names %in% drops)
    X <- X[,-drop.indices]; W <- W[,-drop.indices]
    rg <- rg[-drop.indices,-drop.indices]
    names <- names[-drop.indices]
  }

  #remove NAs, extreme values.
  r <- matrixGWASQC(X,W,all_ids, is.sim = args$simulation)
  X <- r$clean_X;  W <- r$clean_W; all_ids <- r$snp.list
  if(length(r$dropped_cols) > 1 || r$dropped_cols != 0)
  {
    names <- names[-(r$dropped_cols)]
  }

  #Consider moving this into GWASQC, limit studies with fewer than N samples?
  if(args$scale_n != "")
  {
    N <- as.matrix(data.table::fread(args$scale_n, check.names = TRUE) %>% dplyr::filter(!row_number() %in% r$drops) %>%
      dplyr::filter(unlist(.[,1]) %in% all_ids) %>% quickSort(.,args))
    if(args$drop_phenos != "")
    {
      drops <- unlist(strsplit(args$drop_phenos, split = ",")[[1]])
      drop.indices <- which(colnames(N) %in% drops)
      N <- N[,-drop.indices];
      print(N)
      dim(N)
    }
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
      N <- rbind(N, new_rows) %>% quickSort(.,args)
      stopifnot(all(all_ids == N[,1]))
    }
    N <- as.matrix(N[,-1] %>% apply(., 2, as.numeric)) %>% orderColumnsByName(query.mat = ., ref=names,force.ref = args$trait_names)

    #if(args$drop_phenos != "") {N <- N[,-c(as.integer(drop.indices))] } #Drop }

   if(length(r$dropped_cols) > 1 | r$dropped_cols != 0)
   {
     N <- N[,-r$dropped_cols]
   }
    print(dim(N))
    print(dim(W))
    W <- W * (sqrt(N)) #I have been doing this wrong the whole time....
  }
  if(args$genomic_correction != "")
  {
    message("Including genomic correction in factorization...")
    GC <- readInLambdaGC(args$genomic_correction,X, names) #Check this
    X <- X * (1/sqrt(GC)) # we don't weight by this, we actually adjust the effect sizes by this.
  }

   message("Reading in covariance structure from sample overlap...")
  covar.dat <- SampleOverlapCovarHandler(args, names, X)
  return(list("X" = X, "W" = W, "ids" = all_ids, "trait_names" = names, "C" = covar.dat$C,
              "W_c" = covar.dat$W_c, "rg"=rg, "C_block"=covar.dat$C_block))

}

#' Wrapper to get all the important data extracted and used for cohort overlap covariance adjjustment
#'
#' @param args
#' @param names
#' @param X
#'
#' @return list containing Whitening matrix W_c inverse, Block versino of the C matrix and unblocked version of C matrix.
#' @export
#'
SampleOverlapCovarHandler <- function(args, names, X)
{
  #Just read in the file
  C <- readInCovariance(args$covar_matrix, names)
  #If we are scaling by the sample standard deviation, assuming we use LDSC input
  if(args$sample_sd != "")
  {
    message("Verify that the name order is correct. Should be...")
    sd.df <- fread(args$sample_sd)
    sd.scaling <- 1/(as.matrix(sd.df$V2) %*% t(as.matrix(sd.df$V2)))
    diag(sd.scaling) <- 1 #make it correlation matrix
    C <- C*sd.scaling
  }
  #perform covariance matrix shrinkage

  #option 1- perform shrinkage with block adjustment
  adjusted.C <- linearShrinkLWSimple(C, args$WLgamma)
  whitening.dat <- buildWhiteningMatrix(adjusted.C, ncol(X),blockify = args$block_covar)
  return(list("W_c" =  whitening.dat$W_c, "C_block"=whitening.dat$C_block,"C" = adjusted.C))
}



#' Select the K to initialize with. A few different options, made for testing.
#'
#' @param args a list with settings; key settings are K desired (either a number or a method)
#' @param X input data effect sizes
#' @param W input data weights
#' @param svs (optional) the singular values
#'
#' @return a number of K to initialize with
#' @export
#'
#' @examples
selectInitK <- function(args,X,W, svs = NULL)
{
  #Optiosn now are: MAX (default), KAISER, GD
  #GD method: from https://arxiv.org/pdf/1305.5870.pdf
  #Kaiser: from ??? TODO
  #MAX: m-1
  #M/2: M/2

  if(is.null(svs) & args$K %in% c("KAISER", "GD"))
  {
    #We've already done spectral decomp, no need to do it again...
    svrun <- svd(X*W)
    svs <- svrun$d

  }

  if(is.numeric(args$K) & args$K !=0)
  {
    message("User has specified a number of factors initially- this will be given preference.")
    message("No selection method will be used.")
    return(args$K)
  } else
  {
    k = switch(
      args$K,
      "MAX"= ncol(X)-1,
      "GD"= PCAtools::chooseGavishDonoho(X*W,  var.explained=svs^2, noise=1), #see https://rdrr.io/github/kevinblighe/PCAtools/man/chooseGavishDonoho.html
      "KAISER"= sum(svs^2 > mean(svs^2)),
      "K/2"=ceiling(ncol(X)/2),
      "K-2"=ceiling(ncol(X)/2)
    )
    message("Proceeding with initialized K of ", k)

  }

  return(k)
}

readInSettings <- function(args)
{
 option <- list()
	#careful with this environment passed by reference- we do object copying versions so need to be consisgtent.
	#larger changes required if you're going to use this-- need to change some of the BIc functions
	# option <- listenv::listenv()
  option$calibrate_k <- args$calibrate_k
  option[['K']] <- args$nfactors
  option[['iter']] <- args$niter
  option[['convF']] <- 0
 option[["nsplits"]] <- as.numeric(args$cores)
  option[['conv0']] <- args$converged_obj_change
  option[['ones']] <- FALSE
  option[["plots"]] <- args$overview_plots
  option[['disp']] <- FALSE
  #F matrix initialization
  option[['f_init']] <- args$init_F
  option[['epsilon']] <- as.numeric(args$epsilon)
  option[['u_init']] <- args$init_L
  option[["preinitialize"]] <- FALSE
  option[['carry_coeffs']] <- FALSE
  option[["glmnet"]] <- FALSE
  option[["parallel"]] <- FALSE
  option[["fastReg"]] <- FALSE
  option[["ridge_L"]] <- FALSE
  option[['debug']] <- FALSE #args$debug
  option[["subsample"]] <- args$subsample
  #option[["gls"]] <- ifelse(args$covar_matrix != "", TRUE, FALSE)
  option$gls <- FALSE
  option[["covar"]] <- args$covar_matrix
  option$auto_grid_search <- args$auto_grid_search
  option$sort <- args$sort
  if(args$simulation)
  {
    message("specifying minimum K with simulation.")
    option$Kmin <- args$nfactors
  }else
  {
    option$Kmin <- 0
  }
  #Experimental
  message("Scaling is off by default")
  #had to turn off for simulations, at least for now...
  option$scale <- FALSE
  option$burn.in <- 0
  option$fix.alt.setting <- NA
  option$swap <- FALSE
  option$bic.var <- args$bic_var
  option$svd_init <- args$svd_init
  option$std_y <- args$std_y
  option$param_conv_criteria <- args$param_conv_criteria

  #Internal use only:
  option$actively_calibrating_sparsity <- FALSE
  #This looks alot like object-oriented programming.
  #You should just have this be a proper R object with all the attributes and data you need....


  if(args$regression_method %in% c("penalized", "glmnet", "OLS"))
  {
	 option[["regression_method"]] = args$regression_method #push this through all initializations.
  }else{
	  message("This method isn't recognized. Try penalized or glmnet")
	  quit()
  }
  option[["posF"]] <- args$posF
  option$out <- args$output
  #option$logu <- args$output
  option$logu <- NULL #for not writing to log fie output.
  option[["MAP_autofit"]] <- as.numeric(args$MAP_autofit)
  option$intercept_ubiq <- FALSE
  option$traitSpecificVar <- FALSE
  option$V <- args$verbosity
  option$calibrate_sparsity <- args$scaled_sparsity
  if(args$cores > 1)
  {
    updateLog(paste("Running in parallel on", args$cores, "cores"))
    option[["parallel"]] <- TRUE
  }
  option[["ncores"]] <- args$cores
  option[["fixed_ubiq"]] <- args$fixed_first
  option$std_coef <- args$std_coef
  return(option)
}



#####
## Setting defaults helpful for running elsewhere

defaultSettings <- function(K=0, init.mat = "V")
{
  args <- UdlerArgs()
  args$niter <- 200
  args$uncertainty <- ""

  args$gwas_effects <- ""
  args$nfactors <- K
  args$scale_n <- ""
  args$output <- "/scratch16/abattle4/ashton/snp_networks/scratch/testing_gwasMF_code/matrix_simulations/RUN"
  opath <- "gwasMF"
  args$simulation <- TRUE
  args$sort <- FALSE #b/c default for sims.
  args$converged_obj_change <- 0.001
  args$std_coef <- FALSE
  args$std_y <- FALSE
  if(init.mat == "U")
  {
    args$init_L <- "std"
  }
  args$fixed_first <- TRUE #trying to see if this help- IT DOENS'T really appear to matter very much.
  args
}

DefaultSeed2Args <- function()
{
  args <- list()
  args$std_coef <- FALSE
  args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/B.tsv"
  args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/SE.tsv"
  args$fixed_first <- TRUE
  args$genomic_correction <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/Lambda_gc.tsv"
  args$nfactors <- 15
  args$calibrate_k <- FALSE
  args$trait_names = ""
  args$niter <- 10
  args$alphas <- ""
  args$lambdas <- ""
  args$autofit <- -1
  args$auto_grid_search <- TRUE
  args$cores <- 1
  args$IRNT <- FALSE
  args$weighting_scheme = "B_SE"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/scratch/testing_gwasMF_code/model_selection/bic_autofit/"
  args$converged_obj_change <- 1
  args$scaled_sparsity <- TRUE
  args$posF <- FALSE
  args$init_F <- "ones_eigenvect"
  args$init_L <- ""
  args$epsilon <- 1e-8
  args$verbosity <- 1
  args$scale_n <- ""
  args$MAP_autofit <- -1
  args$auto_grid_search <- FALSE
  args$regression_method = "penalized"
  args$converged_obj_change <- 0.001
  args$sort <- TRUE
  args
}
YuanSimEasy <- function()
{
  args <- list()
  args$std_coef <- FALSE
  args$covar_matrix = ""
  args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/yuan_simulations/Input_tau100_seed1_X.txt"
  args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/yuan_simulations/Input_tau100_seed1_W.txt"
  args$fixed_first <- TRUE
  args$genomic_correction <- ""
  args$overview_plots <- FALSE
  args$nfactors <- 5
  args$calibrate_k <- FALSE
  args$trait_names = ""
  args$niter <- 100
  args$simulation <- TRUE
  args$alphas <- ""
  args$lambdas <- ""
  args$autofit <- -1
  args$auto_grid_search <- FALSE
  args$cores <- 1
  args$IRNT <- FALSE
  args$weighting_scheme = "B_SE"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/yuan_simulations/"
  args$converged_obj_change <- 1
  args$scaled_sparsity <- TRUE
  args$posF <- FALSE
  args$init_F <- "ones_eigenvect"
  args$init_L <- ""
  args$epsilon <- 1e-8
  args$verbosity <- 1
  args$sort <- FALSE
  args$scale_n <- ""
  args$MAP_autofit <- -1
  args$auto_grid_search <- FALSE
  args$regression_method = "penalized"
  args$converged_obj_change <- 0.05 #this is the percent change from one to the next.
  args$prefix <- ""
  args$bic_var <- "mle"
  args$svd_init <- FALSE
  args
}
UdlerArgs <- function()
{
  args <- list()
  args$sort <- TRUE
  args$std_coef <- FALSE
  args$covar_matrix = ""
  args$gwas_effects <-"/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/beta_signed_matrix.tsv"
  args$uncertainty <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/se_matrix.tsv"
  args$fixed_first <- TRUE
  args$genomic_correction <- ""
  args$overview_plots <- FALSE
  args$nfactors <- 57
  args$calibrate_k <- FALSE
  args$trait_names = ""
  args$niter <- 150
  args$alphas <- ""
  args$lambdas <- ""
  args$autofit <- -1
  args$auto_grid_search <- TRUE
  args$cores <- 1
  args$svd_init <- TRUE
  args$IRNT <- FALSE
  args$weighting_scheme = "B_SE"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/udler_original/bic_version/"
  args$converged_obj_change <- 1
  args$scaled_sparsity <- TRUE
  args$posF <- FALSE
  args$std_y <- FALSE
  args$init_F <- "ones_eigenvect"
  args$init_L <- ""
  args$epsilon <- 1e-8
  args$verbosity <- 1
  args$scale_n <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/sample_counts_matrix.tsv"
  args$MAP_autofit <- -1
  args$auto_grid_search <- FALSE
  args$regression_method = "glmnet"
  args$converged_obj_change <- 0.05 #this is the percent change from one to the next.
  args$prefix <- ""
  args$bic_var <- "mle"
  args$param_conv_criteria <- "BIC.change"
  args
}

#Try another selection method..
#https://rdrr.io/github/kevinblighe/PCAtools/src/R/randomMethods.R
#kevinblighe/PCAtools
#I take no credit for this work whatsoever.

#' Filling in settings, all relics of a past run.
#' TODO: basically get rid of all this
#'
#' @param curr.args
#'
#' @return updated arguments
#' @export
#'
fillDefaultSettings <- function(curr.args)
{
  curr.args$std_coef <- FALSE
  curr.args$sort <- TRUE #what is this doing? sorts the snps I think...
  curr.args$calibrate_k <- FALSE
  curr.args$alphas <- ""
  curr.args$lambdas <- ""
  curr.args$autofit <- -1
  curr.args$auto_grid_search <- TRUE
  curr.args$cores <- 1
  curr.args$svd_init <- TRUE
  curr.args$IRNT <- FALSE
  curr.args$scaled_sparsity <- TRUE
  curr.args$posF <- FALSE
  curr.args$std_y <- FALSE
  curr.args$init_F <- "ones_eigenvect"
  curr.args$init_L <- ""
  curr.args$MAP_autofit <- -1
  curr.args$auto_grid_search <- FALSE
  curr.args$regression_method = "glmnet"
  curr.args$prefix <- ""
  curr.args$scale_n <- ""
  curr.args$output <- curr.args$outdir
  curr.args$simulation <- FALSE
  curr.args$std_coef <- FALSE
  curr.args$std_y <- FALSE
  #TODO: add some checks for missing or cnflicing parameters
  curr.args
}


#' Quick readout to report settings, useful for debugging
#'
#' @param argsin the current arguments in.
#'
#' @return
#' @export
#'
#' @examples
writeRunReport <- function(argsin)
{
  message("Running GLEANER, with settings as follows:")
  message("")
  message("--------------- INPUT FILES ---------------")
  message("Effect sizes: ", basename(argsin$gwas_effects))
  message("Uncertainty estimates: ", basename(argsin$uncertainty))
  message("Cohort overlap adjustment: ", basename(argsin$covar_matrix))
  message("       Block distance: ", as.character(argsin$block_covar))
  message("       Shrinkage factor: ", as.character(argsin$WLgamma))

  message("Genomic correction terms: ", basename(argsin$genomic_correction))
  message("Z-score sample standard deviation: ", basename(argsin$sample_sd))
  message("")

  message("--------------- INPUT SETTINGS ---------------")
  message("BIC convergence criteria: ", argsin$param_conv_criteria)
  message("BIC method: ", argsin$bic_var)
  message("K init: ", argsin$nfactors)
}

#Projection helper functions
#These come from Yuan He's code, and have been modified for my use here....

#Weight Some matrix by W. Flexible to take either F or X
weightEffects <- function(X,W, z_score = FALSE)
{
  #We presume if W is Ses, that W <- 1/SE
  #X cam 
  #Return item as a list
  if(is.null(W))
  {
    message("WARNING: null W")
    return(lapply(1:nrow(X), function(i) X[i,]))
           
  } else if(all(dim(X) == dim(W)) | z_score) {
    #return(lapply(1:nrow(X), function(i) X[i,]*W[i,]))
    return(split(X*W,seq(nrow(X))))
    
  } else { 
    N <- nrow(W)
    #Code nippet from YUAn
    weighted_data = list()
      for(n in seq(1,N)){
        weighted_data[[n]] = list()
        w = as.numeric(W[n, ])
        weighted_data[[n]] = diag(w) %*% FactorM_optimal
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
  if(is.na(trait_names))
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
  compositions = apply(FactorM_optimal, 2, function(x) traits[which(x>thresh)])
  
  return(list("F"=FactorM_optimal, "composition"=compositions, "trait"=unlist(traits)))
}

#Now actually do the projection
lmProjection <- function(factors, projection.target, W, trait.list, f.composition)
{
  #toy W:
  if(is.null(W))
  {
    W <- matrix(rep(1, nrow(projection.target)* ncol(projection.target)), 
                nrow = nrow(projection.target), ncol = ncol((projection.target)))
  }
  message("Weighting F")
  weighted.F <- weightEffects(factors, W)
  message("Weighting X")
  weighted.X <- weightEffects(projection.target, W)
  message("Prepping object")
  process.obj <- lapply(1:length(weighted.X), function(i) list("X" = weighted.X[[i]], "F" = weighted.F[[i]]))
  return(run_regression(process.obj, trait.list, f.composition))
}

#Taken directly from YUAN
run_regression <- function(weighted_data, tissues, compositions){
  options(warn=-1)
  message("performing regression...")
  result = lapply(weighted_data, function(x) compute_p(x, tissues, compositions))
  pValues = c()
  Betas   = c()
  SE = c()
  for(r in result){
    Betas = rbind(Betas, r[[1]])
    pValues = rbind(pValues, r[[3]])
    SE = rbind(SE, r[[2]])
  }
  
  Betas = as.data.frame(Betas)
  pValues = as.data.frame(pValues)
  SE = as.data.frame(SE)
  return(list("B"=Betas,"SE" = SE, "P" =pValues))
}

#Taken, with some tweaks

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
    factor_exist = which(sapply(compositions, function(x) length(intersect(data_exist, x))>length(x)/2))
    
    dataPoint[[1]][invalid_idx] = 0
    if(length(dataPoint[[1]]) != length(tissues)){
      dataPoint[[1]] = t(dataPoint[[1]])
    }
    dataPoint[[2]] = dataPoint[[2]][, factor_exist]
    if(length(dataPoint[[2]]) == length(tissues)){
      dataPoint[[2]] = t(t(dataPoint[[2]]))
    }
    dataPoint[[2]][invalid_idx,] = 0
    
    # remove the shared factor if there are other factors become co-linear with it (because of NA values)/other factors that are also fully loaded on every trait.
    #modified so its not strictly 0, but some non-zero elements allowed (to value of 1e-5)
    if((1 %in% factor_exist) & (length(factor_exist) > 1)){
      temp_d2 = as.matrix(dataPoint[[2]][, 2:length(factor_exist)])
      if(sum(abs(dataPoint[[2]][, 1]) > thresh) == max(apply(temp_d2, 2, function(x) sum(abs(x) > thresh)))){
        dataPoint[[2]] = dataPoint[[2]][, seq(2, length(factor_exist))]
        factor_exist = factor_exist[factor_exist!=1]
      }
    }
    
    if(length(factor_exist) == 0){
      beta = rep(0, length(compositions))
      pv = rep(-1, length(compositions))
      return(list(beta, pv))
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
  if(length(pv) < length(compositions)){
    beta_formatted = rep(0, length(compositions))
    pv_formatted = rep(-1, length(compositions))
    beta_formatted[factor_exist] = beta
    pv_formatted[factor_exist] = pv
    beta = beta_formatted
    pv = pv_formatted
  }
  return(list("B"=as.numeric(beta), "SE"=as.numeric(se), "P"=as.numeric(pv)))
}


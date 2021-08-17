# Functions to assess the performance of benchmarked methods
##Comparison methods:

#Frobenius norm difference from real matrix
frobReconstruction <- function(pred, true)
{
  return(norm(true - pred, "F"))
}

#Helper function to find the optimal correlation permutation of the columns....
#I might just go with Yuan on this one, this is too much... also it doesn't guarantee just one, 
#
factorMapping <- function(tcc_output, nfact_hat, nfact_true)
{
  t <- data.frame(tcc_output) %>% mutate("lhat" = 1:nfact_hat)
  nl <- 1:nfact_true
  colnames(t) <- c(nl, "lhat")
  c <- pivot_longer(t, nl) %>% mutate(val2 = value^2) %>% arrange(-val2, decreasing = TRUE)
  assigned_loading <- c()
  true_loading <- c()
  corr_vals <- c()
  i = 1
  while(length(assigned_loading) < nrow(tcc_output) & i < nrow(c))
  {
    curr = c[i,]
    if(!(curr$lhat %in% assigned_loading) & !(curr$name %in% true_loading)) #if you haven't already matched a factor
    {
      assigned_loading <- c(assigned_loading, as.character(curr$lhat))
      true_loading <- c(true_loading, curr$name)
      corr_vals <- c(corr_vals, curr$value)
    }
    i = i + 1
  }
  missing <- nl[!(nl %in% true_loading)]
  a <- cbind(rep("-",length(missing)), missing, rep(0,length(missing)))
  final <- cbind(assigned_loading, true_loading, corr_vals)
  t <- data.frame(rbind(final, a))
  t$corr_vals <- as.numeric(as.character(t$corr_vals))
  return(t)
}

#function to match up columns by their correlation. Probably just bad janky code.
#Get teh correlations and pull them out
topCorCols <- function(pred, true)
{
  
  m_cor <- cor(pred, true)
  return(factorMapping(m_cor, ncol(pred), ncol(true)))
}

getAverageCorrelation <- function(pred, true, omit_zeros = FALSE)
{
  t<- topCorCols(pred, true)
  if(omit_zeros)
  {
    r <- t$corr_vals[t$corr_vals != 0]
    return(mean(abs(r)))
  }
  mean(abs(t$corr_vals))
}


############
##### From Yuan
getCoordinatedCorrelation <- function(fp, lp, trueF, trueL, rankK){
  #Fill in missing column
  if(ncol(fp) < ncol(trueF)){
    dif = ncol(trueF) - ncol(fp)
    fp = cbind(fp, matrix(rep(0, nrow(fp) * dif), ncol = dif))
    lp = cbind(lp, matrix(rep(0, nrow(lp) * dif), ncol = dif))
    sig_hits = cbind(sig_hits, matrix(rep(0, nrow(sig_hits) * dif), ncol = dif))
  }
  rankK = ncol(trueF)
  suppressWarnings(library('combinat'))
  ordering = permn(seq(1,rankK))
  #get the correlation matrix in every order possible
  f_cor = rep(0, length(ordering))
  for(ord in seq(1, length(ordering))){
    f_cor[ord] = cor(as.vector(trueF), as.vector(fp[,ordering[[ord]]]))
  }
  
  l_cor = rep(0, length(ordering))
  for(ord in seq(1, length(ordering))){
    l_cor[ord] = cor(as.vector(trueL), as.vector(lp[,ordering[[ord]]]))
  }
  
  #select which one gives the greatest overall correlation across both f and l
  ord_sum = f_cor + l_cor
  ord = which.max(ord_sum)
  lp = lp[,ordering[[ord]]]
  fp = fp[,ordering[[ord]]]
  sig_hits = sig_hits[, ordering[[ord]]]
  
  lp = lp / matrix(rep(apply(lp, 2, function(x) max(abs(x))), nrow(lp)), nrow = nrow(lp), byrow = T)
  fp = fp / matrix(rep(apply(fp, 2, function(x) max(abs(x))), nrow(fp)), nrow = nrow(fp), byrow = T)
  colnames(fp) = seq(1,ncol(fp))
  rownames(fp) = seq(1, nrow(fp))
  
  fp[is.na(fp)] = 0 
  lp[is.na(lp)] = 0
  
  ret <- list()
  ret$l_corr = cor(as.vector(trueL), as.vector(lp))
  ret$f_corr = cor(as.vector(trueF), as.vector(fp))
  return(ret)
}
#####



#This function takes a list of (F, L) tuples and calculates performance across them all
assessMultiple <- function(results, trueF, trueL, names)
{
  results 
  for(i in 1:length(results))
  {
    
  }
}

#This function takes an F, L and calculates its performance
#TODO: perhaps we don't want average correlation, but correlation of each factor?
assessFactorization <- function(f_hat, l_hat, x_hat, F, L)
{
  ret <- list()
  corrs <- getCoordinatedCorrelation(f_hat, l_hat, F, L, ncol(L))
  #ret$l_corr <- getAverageCorrelation(l_hat, L) 
  #ret$f_corr <- getAverageCorrelation(f_hat, F)
  ret$l_corr <- corrs$l_corr
  ret$f_corr <- corrs$f_corr
  ret$frob <- frobReconstruction(x_hat, L %*% t(F))
  return(ret)
  
}
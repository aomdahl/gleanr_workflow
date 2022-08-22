copheneticAssessment <- function(M){
  #compute the consensus matrix C, allowing samples to be assigned to one factor
  C = list();
  D = nrow(M[[1]])
  
  for(m in 1:length(M)){
    assignment = apply(M[[m]], 1, which.max)
    ag_rep_rows = matrix(rep(assignment, D), ncol=D)
    C[[m]] = (ag_rep_rows == t(ag_rep_rows)) * 1;
  }
  #take the average to get the fraction of consistent assignment
  Cbar = aaply(laply(.data=C, .fun=as.matrix), c(2, 3), mean)
  
  #compute the correlation between original distances and the cophenetic distances from a hierachical clustering based on average linkage
  if(sum(Cbar != 1) == 0){
    return(1)
  }
  Ds = Cbar - diag(rep(1,D));
  if(sum(Ds) == 0){
    return(1)
  }
  d1 = dist(Ds)
  hc = hclust(d1, "ave")
  d2 = cophenetic(hc)
  coph = cor(d1, d2)
  
  return(coph)
}

correlationAssessment <- function(M)
{
  
  factor_corrs = lapply(M, cor)
  corr_frobs = lapply(factor_corrs, function(x) norm(x, 'F'))
  return(mean(unlist(corr_frobs)))
  
}
scaledCorrelationAssessment <- function(M) {
  factor_corrs = lapply(M, cor)
  print(ncol(M[[1]]))
  corr_frobs = lapply(factor_corrs, function(x) norm(x, 'F'))
  return(mean(unlist(corr_frobs))/ncol(M[[1]]))
}

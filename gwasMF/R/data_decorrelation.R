######
## Functions for adjusting for covariance structure by building a decorrelation matrix (to be used as in GLS)
######

#' Force a covariance matrix to be block diagonal
#' This matches our assumption that studies share covariance with a block diagonal structure (i.e. limited numbers are grouped)
#'
#' @param covar The actual covariance matrix; a matrix object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
blockifyCovarianceMatrix <- function(covar,...)
{
  blocks <- create_blocks(covar,...)
  rownames(covar) <- colnames(covar)
  blockify(covar, blocks)
}

#From Guanghao fast asset
create_blocks <- function(cormat, cor_thr=0.2){
  # Hierarchical clustering
  corrdist <- as.dist(1-abs(cormat))
  hc <- hclust(corrdist)
  htree <- cutree(hc, h=1-cor_thr) # Criterion: corr>0.2
  block <- as.integer(names(table(htree))[table(htree)>=2])
  block <- lapply(block, function(x) names(htree)[htree==x])

  return(block)
}

blockify <- function(cormat, blocks)
{
  ret.mat <- diag(diag(cormat))
  colnames(ret.mat) <- colnames(cormat); rownames(ret.mat) <- rownames(cormat);
  if(length(blocks) != 0)
  {

    for(b in blocks)
    {
      if(length(b) >= 2)
      {
        for(i in b)
        {
          for(j in b)
          {
            ret.mat[i,j] <- cormat[i,j]
          }
        }
      }
    }
  }
  return(ret.mat)

}


#buildWhiteningMatrix(decorrelate)
buildWhiteningMatrix<- function(covar, dim, blockify = FALSE,...)
{
  if(is.null(covar)){
    message("buildWhiteningMatrix- decorrelating with an identity matrix.")
    return(diag((dim))) #Just return an identity matrix, we aren't going to do anything....
  }
  if(blockify)
  {
    covar <- blockifyCovarianceMatrix(covar,...)
  }
  if(!isSymmetric(covar, tol=1e-3))
  {
    covar <- as.matrix(Matrix::forceSymmetric(covar))
  }
  #make sure covar is symmetric and pd.
  if(!lqmm::is.positive.definite(covar, tol = 1e-3))
  {
    message("C is not PD, adjusting now...")
    covar <- lqmm::make.positive.definite(covar)
  }

  #I think this was wrong...
  #solve(t(chol(covar)))
  solve((chol(covar)))
}

#Account for whitening, if sppecified:
#Covar- the covariance you want to account for
#If this is null, no change made
#whitener- don't calculate the whitening transofrm new, assume its been done already.
#adjustMatrixAsNeeded(diag(w) %*% X, decorrelate, whitener = precalc.U)
adjustMatrixAsNeeded <- function(X, covar, whitener = NULL)
{
  X <- as.matrix(X)
  if(!is.null(covar))
  {
    if(is.null(whitener))
    {
      U_t_inv <- buildWhiteningMatrix(covar, ncol(X))
    }else
    {
      U_t_inv <- whitener
    }
    #check the dimensions. Only need the double transpose if the dimensions aren't aligned.
    #In the case that they are, just leave alone
    if(nrow(X) == nrow(U_t_inv))
    {
      Xin <- U_t_inv %*% (X)
    }else #this is the joined matrix case, some transposition here to make it work.
    {
      message("NEeding to transposes")
      print(head(U_t_inv))
      print(head(X))
      Xin <- t(U_t_inv %*% t(X))
      #Sanity check
      C2 <- var(X)
      U_t_inv <- solve(t(chol(C2)))
      Xin2 <- t(U_t_inv %*% t(X))
      if(any(round(var(Xin2), digits = 3) != diag(nrow(covar))))
      {
        message("Possible dangerous territory, error with whitening....")
        non.compliant <- which(round(var(Xin2), digits = 3) != 1 & round(var(Xin2), digits = 3) != 0)
        print(var(Xin2)[non.compliant])
      }
    }
    stopifnot(all(dim(Xin) == dim(X)))

    return(Xin)
  } else
  {
    if(!is.null(whitener))
    {
      message("No covariance structure passed, but whitener is specified.")
      message("NO CORRECTION BEING MADE.")
    }
    return(X)
  }
}

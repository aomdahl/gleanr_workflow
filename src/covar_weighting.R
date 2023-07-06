#U must be a list of all the choleskyt decomposition matrices we are going to operate on.
#This is helpful so we only do this once.
#@param U.inv.l a list of all the whitening matries
preGLSScaleBenchmark <- function(X.mat, U.inv.l)
{
  library(Matrix)
  
  #2 ways to do this, either as an lapply typ eof thing or a matrix product. Should benchmark both and see which is faster
  #Here, we are correcting for one source of covariance at a time (either cohort overlap or LD)
  X.mat <- rbind(X,X,X,X,X,X,X,X,X,X,X,X)
  #Check dimensions
  
  U.inv.l <- lapply(1:nrow(X.mat), function(i) {x <- diag(X.mat[i,]); diag(x) <- rnorm(10); x})
  #Method 1: Lapply
  #For each row in X, adjust it by the whitening matrix given by U.inv.l. We assume that these correspond by row, but need a way to enforce this.
  simple <- system.time({adjusted.x <- do.call("rbind", lapply(1:nrow(X.mat), 
                                                     function(i) t(adjustMatrixAsNeeded(t(X.mat[i,]), 
                                                                                      U.inv.l[[i]], whitener = U.inv.l[[i]]))))})
  
  #Method 1:
  #Expand U to a sparse diagonal vector
  #Expand X to a long vector, where rows are stacked (not columns)
  complex <- system.time({
    x.v <- as.vector(as.matrix(t(X.mat)))
    o <- Matrix::bdiag(U.inv.l)
    #Matrix::kronecker()
    product <- o %*% x.v
    adjusted.x.mat <- as.matrix(matrix(product, nrow = nrow(X.mat), byrow = TRUE))
  })
  stopifnot(all(adjusted.x == adjusted.x.mat))
  print(simple)
  print(complex)
  message("Imn conclusion, it appears that the matrix version runs about 2x as fast. ")
}


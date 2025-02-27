#### Function to evaluate error between predicted L and F and true L and F
#### Taken from Yuan He's code with some modification.
unitNorm <- function(x)
{
  x/norm(x, "2")
}
unitNorms <- function(M)
{
  ret <- apply(M, 2, function(x) unitNorm(x))
  ret[is.na(ret)] <- 0
  ret
}

#Changes:
  #Sig_hits not considered
#lp <- pred.loadings
#fp <- pred.factors
#
#evaluteFactorConstruction(true.loadings, true.factors, pred.loadings, pred.factors)
if(FALSE)
{
  trueL=true.loadings
  trueF = true.factors
  lp= pred.loadings
  fp=pred.factors
}
#as.matrix(true.v), as.matrix(FactorM)
#evaluateSingleMatrixConstruction(as.matrix(true.v), as.matrix(ret$V))
evaluateSingleMatrixConstruction <- function(truem, predm)
{
  truem <- as.matrix(truem); predm <- as.matrix(predm);

  #fill in missing columns with 0s
  if(ncol(predm) < ncol(truem)){
    #message("missing factors, will set to 0...")
    dif = ncol(truem) - ncol(predm)
    predm = cbind(predm, matrix(rep(0, nrow(predm) * dif), ncol = dif))
  }
  #What if the new one learns
  if(ncol(truem) < ncol(predm)){
    #message("missing factors, will set to 0...")
    dif = ncol(predm) - ncol(truem)
    truem = cbind(truem, matrix(rep(0, nrow(predm) * dif), ncol = dif))
  }


  #fill in NAs
  predm[is.na(predm)] = 0

  K = ncol(truem) #previously was truem, had to change because
  suppressWarnings(library('combinat'))
  if(K > 9)
  {
    message("There are too many latent factors for comprehensive examination of options. Do a greedy search.")
    return()
  } else
  {
    ordering = permn(K)
    sign.ordering <- gtools::permutations(2, K, c(1,-1), repeats.allowed = TRUE)
    all.l <- getallR2(sign.ordering, ordering, predm, truem)
    return(all.l[which.max(all.l)])
  }
}

#updated.mat <-fillWithZeros(trueF, fp, lp=lp)
fillWithZeros <- function(trueF, fp, lp = NULL)
{
  trueF <- as.matrix(trueF)
  fp <- as.matrix(fp)
  if(ncol(fp) < ncol(trueF)){
    #message("missing factors, will set to 0...")
    dif = ncol(trueF) - ncol(fp)
    fp = cbind(fp, matrix(rep(0, nrow(fp) * dif), ncol = dif))
    if(!is.null(lp))
    {
      lp = cbind(lp, matrix(rep(0, nrow(lp) * dif), ncol = dif))
    }
    
  }
  #Another special case I found....
  if(nrow(fp) == 0)
  {
    fp <- matrix(0, nrow = nrow(trueF), ncol = ncol(trueF))
  }
  return(list("fp" = fp, "lp" = lp))
}

#Ordering may be specified as a list of lists, with possible orders
#sel.u, sel.v, res$U, res$V
#evaluteFactorConstruction(true.loadings, true.factors, pred.loadings, pred.factors,unit.scale = FALSE)
#This version requires that L and F have the same ordering.
#May need to circle back on this...
#Its possible
evaluteFactorConstruction <- function(trueL, trueF, lp, fp, easy=TRUE, verbose = FALSE, ordering = NULL,corr.type = "pearson",null.return = FALSE,...){
  if(null.return)
  {
    return(c(NA, NA))
  }
  if(is.null(lp) | length(lp) == 0)
  {
    message("empty table")
    lp <- matrix(0, nrow= nrow(trueL), ncol = nrow(trueL))
  }
  if(is.null(fp) | length(fp) == 0)
  {
    message("empty table")
    fp <- matrix(0, nrow= nrow(trueF), ncol = nrow(trueF))
  }

  #make sure all are matrices

  trueL <- as.matrix(trueL); trueF <- as.matrix(trueF); lp <- as.matrix(lp); fp <- as.matrix(fp)


  #fill in missing columns with 0s
  updated.mat <-fillWithZeros(trueF, fp, lp=lp)
  fp <- updated.mat$fp; lp <- updated.mat$lp

  #fill in NAs
  fp[is.na(fp)] = 0
  lp[is.na(lp)] = 0

  K = ncol(trueF)
  suppressWarnings(library('combinat'))
  if(K > 9 & is.null(ordering))
  {
    message("There are too many latent factors for comprehensive examination of options. Reverting to a greedy search, in which U and V are treated separately.")
    return(c(greedyMaxCorr(trueL, lp,cor.type = corr.type), greedyMaxCorr(trueF, fp,cor.type = corr.type)))
  } else if(easy)
  {
    if(is.null(ordering))
    {
      ordering = permn(K)
    }

    sign.ordering <- gtools::permutations(2, K, c(1,-1), repeats.allowed = TRUE)
    all.l <- getallR2(sign.ordering, ordering, lp, trueL,...)
    all.f <- getallR2(sign.ordering, ordering, fp, trueF,...)
    tot = all.l + all.f
    if(all(is.na(all.l)) | all(is.na(all.f)))
    {
      return(c(0,0))
    }
    return(c(all.l[which.max(tot)],all.f[which.max(tot)]))

  }
  else {
    message("Doing this the easy way...")
    #get all possible permutations in F and L
    #need to account for differences in possible sign assignments too....
    ordering = permn(K)
    f.perf <- tableCor(fp, trueF, ordering)
    l.perf <- tableCor(lp, trueL, ordering)
    #Look at the sum of correlations overall (taking absolute value)
    #Note that there may be an exceptional case, where they have opposite directions of correlation, and this is missed
    ord_sum = (f.perf$r2 + l.perf$r2)
    ord = which.max(ord_sum)  #which combination is highest overall
    #Ensure the sign choice is the same - if one is a switched sign, the other must be too....
    if(!(all(f.perf$signs[[ord]] == l.perf$signs[[ord]])) & !(all(f.perf$signs[[ord]] == -1*l.perf$signs[[ord]]))) #also acceptable if they are just inverses of each other, that's an easy fix.
    {
      message("Sign assignments aren't consistent")
      message("Comparing an ordered approach...")
      #Here, fit L (done above), then use those orders on F
      start.pref <- list("l.first" =  tableCor(fp, trueF, ordering, sign.prior = l.perf$signs),
                                "f.first" =  tableCor(lp, trueL, ordering, sign.prior = f.perf$signs))
      #Determine which of these yields the best
      max.choices = (list("l.first" = max(l.perf$r2 + start.pref$l.first$r2), "f.first" = max(f.perf$r2 + start.pref$f.first$r2)))
      #1 is lfirst, 2 is f.first
      #Choose the maximizing one.
      max.fit = start.pref[[which.max(max.choices)]]
      #get the order index of the maximizing one
      ord.second <- which.max(max.fit$r2)
      if(names(which.max(max.choices)) == "l.first")
      {
        l_corr = l.perf$r2[[ord.second]]
        f_corr = start.pref$l.first$r2[[ord.second]]
      } else if (names(which.max(max.choices)) == "f.first"){
        f_corr = f.perf$r2[[ord.second]]
         l_corr = start.pref$f.first$r2[[ord.second]]
      }
      else
      {
        message("ERROR")
        l_corr = NA
        f_corr = NA
      }

    } else {
      l_corr = l.perf$r2[[ord]]
      f_corr = f.perf$r2[[ord]]
    }
    return(c(l_corr, f_corr))
  }

}

#note- there are 2 approaches here:
#1 we try each sign combo for each max combo. O(k^2*k!), very expensive for large K, but cheaper for small
#2 we first identify the best sign for each, then get the sums at the end (O(k^3 _ k!), still expensive but much cheaper

#function
bestR2BySign <- function(pred.v, true) #takes in a full vector version.
{
  library(RcppAlgos)
  REF <- c(1,-1)
  K <- ncol(true)
  N <- nrow(true)
  true.v <- as.vector(true)
  #remove the last one, its all negatives
  sign.options <- RcppAlgos::permuteGeneral(REF, K, TRUE)
  r.choices <- c()
  for(i in 1:nrow(sign.options))
  {
    sign.expand <- unlist(lapply(sign.options[i,], function(x) rep(x, N)))
    r.choices <- c(r.choices, r2Signed(pred.v, sign.expand, true.v))
  }
  best <- which.max(r.choices)
  return(list("r2" = max(r.choices), "signs"= sign.options[best,]))
}

r2Signed <- function(pred.v, signs, true.v)
{
  return(cor(pred.v * signs, as.vector(true.v))^2)
}


#wrapper function to check all the possible orderings.
tableCor <- function(pred, true, ordering, sign.prior = NULL){
  f_cor = rep(0, length(ordering))
  f_signs = list()
  for(ord in seq(1, length(ordering))){
    curr.order <- as.vector(pred[,ordering[[ord]]]) #for the current ordering, which combination of signs yields the best fit?
    if(is.null(sign.prior))
    {
      best.sign <- bestR2BySign(curr.order, true)
    }else
    {
      best.sign <- list("r2"=r2Signed(curr.order, sign.prior[[ord]], true), "signs"=sign.prior[[ord]])
    }

    f_cor[ord] = best.sign$r2
    f_signs[[ord]]=best.sign$signs
  }
  return(list("r2" = f_cor, "signs"=f_signs))
}

#getallR2(sign.ordering, ordering, predm, truem)
getallR2 <- function(sign.ordering, ordering, pred, actual, unit.scale = FALSE)
{
  nreps = nrow(actual)
  all.options.factors <- apply(sign.ordering, 1, function(s) 
    {smat <- do.call("rbind", lapply(1:nreps, function(x) s)); lapply(ordering, function(x) pred[,x] * smat)})
  res.mat <- matrix(NA,nrow(sign.ordering),length(ordering))
  #For each of those, estimate the R2
  if(unit.scale){actual <- unitNorms(actual) }
  for(signopt in 1:nrow(sign.ordering))
  {
    if(unit.scale)
    {
      
      #res.mat[signopt,] <- unlist(lapply(all.options.factors[[signopt]], function(x) cor(as.vector(unitNorms(x)), as.vector(actual))^2))
      #Updating to handle 0 cases and not throw errors
      res.mat[signopt,] <- unlist(lapply(all.options.factors[[signopt]], function(x) customCorr(as.vector(unitNorms(x)), as.vector(actual))^2))
       
    }else
    {
      #res.mat[signopt,] <- unlist(lapply(all.options.factors[[signopt]], function(x) cor(as.vector(x), as.vector(actual))^2))
      #Updating to handle 0 cases and not throw errors
      res.mat[signopt,] <- unlist(lapply(all.options.factors[[signopt]], function(x) customCorr(as.vector(x), as.vector(actual))^2))
    }
    
  }
  return(res.mat)
}


#####5/31
###A greedy implementation to evaluate factor construction, one matrix at a time.
#The main recursive call- requires
#Correlation matrix between remaining factors
#R2 matrix of remaining factors
#signs- vector of the signs assigned to each column
#pairs-a list where each entry is a tuple of (reference vect, pred vect) that match to each other
#returns signs and pairs at the end.
#recurseGreedyMap(r.mat,r2, c(dir), list(dt))
recurseGreedyMap <- function(cor, r2, signs, pairs)
{
  dt <- which(r2 == max(r2), arr.ind = TRUE)[1,] #just go with the first one...
  dir <- sign(cor[dt[1], dt[2]])
  if(all(r2 == 0) | ncol(r2) == 1)
  {
    return(list("signs"=signs, "pairs" = pairs))
  }
  pairs[[length(pairs) + 1]] <- dt
  
  cor[dt[1],] <- 0; cor[,dt[2]] <- 0
  r2[dt[1],] <- 0; r2[,dt[2]] <- 0
  recurseGreedyMap(cor,r2, c(signs, dir), pairs)
}



#' Helper to ensure correlation is calculated safely
#'
#' @param x either a vector or matrix of values to check, usually the prediction/estimate
#' @param y either a vector or matrix of values to check, usually the reference
#'
#' @return a list with items PASS <logical>, indicating if it passes, and ret_cor <double> specifying the correlation to return if PASS is false.
#' @export
#'
#' @examples
corrSafteyChecks <- function(x,y, cor.type = "pearson")
{
  ret_l <- list("PASS"=TRUE, "ret_cor"=NA)
  if(all(is.na(x)))
  {
      #if all entries are NA, set them to 0
      x[is.na(x)] <- 0   
  }
  if(all(is.na(y)))
     {
       #if all entries are NA, set them to 0
       y[is.na(y)] <- 0   
  }
  
  if(all(x == y))
  {
    message("Perfect match!")
    ret_l <- list("PASS"=FALSE, "ret_cor"=1)
  }
  if(is.matrix(x) | is.matrix(y))
  {
    var.x <- apply(x,2,var)
    var.y <- apply(y,2,var)
    if(all(var.x == 0) | all(var.y == 0))
    {
      ret_l <- list("PASS"=FALSE, "ret_cor"=diag(0,ncol(y)))
    }
    if(any(var.x == 0) | any(var.y == 0))
    {
      if(any(var.x == 0) & any(var.y == 0))
      {
        warning("Highly unusual case, please revisit")
        message("EvaluateSimR2.R- case where both have empty colums.")
      }
      #Special case- get the correlation, and make the NA columns all 0
      suppressWarnings(ret.mat <- cor(x,y, method=cor.type))
      ret.mat[is.na(ret.mat)] <- 0
      ret_l <- list("PASS"=FALSE, "ret_cor"=ret.mat)
    }
  }else
  {
    if(all(x == 0) | all(y == 0))
    {
      ret_l <- list("PASS"=FALSE, "ret_cor"=0)
    }
  }
  return(ret_l)
}

#' Correlation between matrix columns optimized for my application
#'
#' @param x first matrix 
#' @param y second matrix
#' @param cor.type desired type (pearson, kendall currentlty accepted)
#'
#' @return matrix of correlation estimates (no p-values)
#' @export
customCorr <- function(x,y, cor.type = "pearson")
{
  #Capturing some extreme edge cases
  is.matrix = TRUE
  saftey.vals <- corrSafteyChecks(x,y, cor.type = cor.type)
  if(is.null(dim(x)) & is.null(dim(y)))
  {
    is.matrix = FALSE
  }
  if(!saftey.vals$PASS)
  {
    return(saftey.vals$ret_cor)
  }else
  {
    if(cor.type == "kendall" & is.matrix)
    {
      r.mat <- t(apply(x, 2, function(a) apply(y, 2, function(b) pcaPP::cor.fk(a,b))))
    }else
    {
      r.mat <- cor(x,y, method=cor.type);
    }
    return(r.mat)
  }


}
#Only does it for one vector, the one you specify first
#greedyMaxCorr(true.v, pred.f)
greedyMaxCorr <- function(true.v, fp, verbose = FALSE, cor.type = "pearson")
{
  ret <- fillWithZeros(true.v, fp)
  r.mat <- customCorr(true.v, ret$fp, cor.type=cor.type);

  #fast way to do this with cor.fk
  
  r.mat[is.na(r.mat)] <- 0
  #r2 <-  cor(true.v, ret$fp,method=cor.type)^2; r2[is.na(r2)] <- 0
  r2 <-  r.mat^2; r2[is.na(r2)] <- 0
  dt <- which(r2 == max(r2), arr.ind = TRUE)[1,]
  stopifnot(max(r2) == r2[matrix(dt,ncol=2)])
  dir <- sign(r.mat[dt[1], dt[2]])
  
  #case- there are multiple maximums. Just pick one to start at.....
  
  #Zero-out the R terms that aren't helping us.....
  r.mat[dt[1],] <- 0; r.mat[,dt[2]] <- 0
  r2[dt[1],] <- 0; r2[,dt[2]] <- 0
  matched.cols <- recurseGreedyMap(r.mat,r2, c(dir), list(dt))

  order.list <- 1:ncol(true.v)
  true.order <- sapply(matched.cols$pairs, function(x) x[1]) %>%  c(., order.list[!(1:ncol(true.v) %in% .)])
  pred.order <- sapply(matched.cols$pairs, function(x) x[2])%>%  c(., order.list[!(1:ncol(true.v) %in% .)])
  pred.signs <-  c(matched.cols$signs, rep(1,(ncol(true.v) - length(matched.cols$signs))))
  #scale
  ##7/02- we don't want to be unit norm scaling here, we save that for the scaled run only.
  #pred.one <- as.matrix(ret$fp[,pred.order] %>% unitNorms(.)) %*% diag(pred.signs)
  #pred.two <- (as.matrix(true.v[,true.order]) %>% unitNorms(.))
  pred.one = matrixSignsProduct(ret$fp[,pred.order], pred.signs)
  pred.two <- (as.matrix(true.v[,true.order]))
  c_ <- stackAndAssessCorr(pred.one,pred.two)
  if(!verbose)
  {
    return(c_)
  }else
  {
    #set the true order to what it was, and the pred order as modified:
    order.f <- order(true.order)
    true.order = true.order[order.f]
    pred.order = pred.order[order.f]
    return(list("corr"=c_, "order.true"=true.order, "order.pred"=pred.order, "signs"=(pred.signs)))
  }
  
}




#' Helper function convert the sign based on a list of signs
#' Needed to catch the special case when signs is length 1
#' @param mat mmatrix to multiply
#' @param signs list of signs to change
#'
#' @return matrrix * diag(signs)
#' @export
#'
#' @examples
matrixSignsProduct <- function(mat, signs)
{
  if(length(signs) == 1)
  {
    #message("Single col case.")
    sc <- as.matrix(mat) * signs
  }else
  {
    sc <- as.matrix(mat) %*% diag(signs)
  }
  sc
}


pseudoProcrustes <- function(A,B, corr.type="pearson")
{
  #ret <- fillWithZeros(A, B)
  #TODO- why do we change the order of both matrices? Dpesn't make sense to me....
  cor.dat <- greedyMaxCorr(A, B, verbose = TRUE, cor.type = corr.type)
  new.A <- as.matrix(A[,cor.dat$order.true]);
  new.B <- matrixSignsProduct(B[,cor.dat$order.pred],cor.dat$signs )
  list("A" = new.A, "B" = new.B, "greedy.match"=cor.dat, "factor.cors" = diag(customCorr(new.A, new.B, cor.type = corr.type)))
}

stackAndAssessCorr <- function(A,B, cor.type = "pearson")
{
  pred.vector <- as.vector(A)
  true.vector <- as.vector(B)
  customCorr(pred.vector,true.vector, cor.type = cor.type)
}

#Taken from https://rpubs.com/mengxu/procrustes_analysis
#procestes.corr <- fullProcrustes(lead, scnd)
procrustes <- function(A, B){
  
  if(ncol(A) ==1 & ncol(B)==1)
  {
    A.normalized <- scale(A)
    B.normalized <- scale(B)
    RSS <- norm(A.normalized - B.normalized,  type = "F")
    return(list(A.normalized = A.normalized, B.normalized = B.normalized, rotation.mtx = 1, B.transformed = B.normalized, RSS = RSS))
  }
  
  # center and normalize A 
  A.centered <- t(scale(t(A), center = TRUE, scale = FALSE))
  A.size <- norm(A.centered, type = "F") / (ncol(A) * nrow(A))
  A.normalized <- A.centered / A.size
  
  # center and normalize B
  B.centered <- t(scale(t(B), center = TRUE, scale = FALSE))
  B.size <- norm(B.centered, type = "F") / (ncol(B) * nrow(B))
  B.normalized <- B.centered / B.size

  # Rotation matrix T 
  #if the matrix is big, just get the first few SVs
  m = B.normalized %*% t(A.normalized)
  if(nrow(m) > 1000 & ncol(m) > 1000)
  {
    message("Large matrix, need to use fast SVD. Please wait...")
    #svd.results <- corpcor::fast.svd(m)
    #Executive decision to just do 20 for speed
    svd.results <- RSpectra::svds(m, 20)
  }else{
    svd.results <- svd(m)
  }
  
  U <- svd.results$u
  V <- svd.results$v
  tr <- V %*% t(U)
  
  # B transformed
  B.transformed <- tr %*% B.normalized
  
  # Error after superimposition
  RSS <- norm(A.normalized - B.transformed,  type = "F")
  
  # Return
  return(list(A.normalized = A.normalized, B.normalized = B.normalized, rotation.mtx = T, B.transformed = B.transformed, RSS = RSS))
}

procrustesVegan <- function(A,B, scale = TRUE)
{
  procrust <- vegan::procrustes(A,B,scale=scale) #Did we want symmetric or not?
  #"symmetric is set to TRUE, both solutions are first scaled to unit variance, giving a more scale-independent and symmetric statistic, often known as Procrustes m2"
  #https://john-quensen.com/tutorials/procrustes-analysis/
  #look at "rotation" and "translation" to get the change
  # Error after superimposition
  #cor.comp <- cor(procrust$Yrot,B) #
  #get.coords.x <- unlist(apply(cor.comp^2, 1, which.max))
  mismatch.potential.sum = max(sum(apply(B,2, function(x) sum(x != 0)) == 0),
                               sum(apply(A,2, function(x) sum(x != 0)) == 0))
  if(mismatch.potential.sum == 1) {mismatch.potential.sum = 2}
  #RSS <- norm(procrust$X - procrust$Yrot,  type = "F")
  RSS <- rrmse(procrust$Yrot,procrust$X) #ORDER IS PRED, TRUE
  #get the signs for verification of order
  order <- apply(procrust$rotation^2, 2, which.max)
  if(any(duplicated(order)))
  {
    dup.val <- order[which(duplicated(order))]
    missing.vals <- (1:length(order))[!(1:length(order) %in% order)]
    dup.indices = which(order %in% dup.val)
    dup.vals = order[dup.indices]
    #So we have too many of those in dup.vals and too few of those in dup.indices
    alt.opts <- unique(c(dup.val, missing.vals))
    option.matrix <- matrix(NA, nrow= length(dup.indices), ncol=length(alt.opts))
    max.opt = -Inf
    max.choice <- NULL
    if(length(dup.indices) > 9)
    {
      message("too many indices mapped to same factor to comprehensively account for this.")
      message("Could try to rectify with existing code, but it should only effect kappa and auprc here, so not going to worry myself about it.")
    }else
    {
      for(possible.choice in combinat::permn(dup.indices))
      {
        choice.matrix <- as.matrix(data.frame(alt.opts,possible.choice))
        curr.opt <- sum((procrust$rotation^2)[choice.matrix])
        if(curr.opt > max.opt) {
          max.opt <- curr.opt
          max.choice <- choice.matrix
        }
      }
      
      #Update the order as appropriate:
      order[max.choice[,2]] <- max.choice[,1]
      
    }
  }

  #check
  cor.comp <- customCorr(procrust$Yrot,B)
  get.coords.x <- unlist(apply(cor.comp^2, 1, which.max))
  mismatch <- sum(get.coords.x != order)
  if(mismatch > mismatch.potential.sum)
  {
    message("Procrustes (single-matrix) mapping to factors isn't clear 1-1. Proceed with caution.")
    #warning("POSSIBLE MAPPING ISSUE with procrustes. DO NOT PROCEED")
  }
  #Is this consistent with what I'd think?
  
  
  signs <- sign(procrust$rotation[matrix(c(1:length(order), order),ncol=2)])

  
  # Return
  return(list(A.normalized = procrust$X, B.normalized = procrust$Yrot, rotation.mtx = T, B.transformed = procrust$Yrot, RSS = RSS,
              "mapping_reorder"=order, "mapping_signs"=signs))
}


# fullProcrustes(lead, scnd)
fullProcrustes <- function(A,B)
{
  #ret <- fillWithZeros(A, B)
  #procrustes(as.matrix(A), as.matrix(ret$fp))
  procrustesVegan(as.matrix(A), as.matrix(B))
}

stackAndAssessNaive <- function(A,B)
{
  #ret <- fillWithZeros(A, B)
  stackAndAssessCorr(A, B)
}

############Looking at Xhats...
rrmse = function(pred,true)
{
  sqrt(sum((pred-true)^2)/sum(true^2))
}
#Could also do this directly in the assessment thing, huh.... urg.
xhatFit <- function(meth, x_hat, x_true, se)
{
  unscaled.methods <- c("PMA2", "backfit_noscale", "flashColumnVar", "factorGo", "GLEANER","GLEANER_glmnet","GLEANER_glmnet_noCovar", "SVD_beta")
  scaled.methods <- c("sSVD", "PCA","SVD", "backfit","ssvd", "sPCA", "PCA_chooseK", "SVD_whiten")
  if(meth %in% unscaled.methods)
  {
    return(list("rrmse"=rrmse(x_hat, x_true), "cor"=stackAndAssessCorr(x_hat, x_true)))
  }else
  {
    #Scale back by the SE
    x_hat_scaled = x_hat * se
    return(list("rrmse"=rrmse(x_hat_scaled, x_true), "cor"=stackAndAssessCorr(x_hat_scaled, x_true)))
  }
  
}

forcedPairedEvaluation <- function(true.u,true.v,pred.u,pred.v)
{
  #procrustes based approach
  #This includes both the U-first approach and my implementation which accounts for both
  if(!(matrixDimsAlign(true.u, pred.u) & matrixDimsAlign(true.v, pred.v))) {
    warning("Matrix dimensions aren't aligned")
    return(NA)
  }
  procrust.dat <- procrustesPairedUV(true.u,true.v,pred.u,pred.v)
  corr.based <- greedyMaxPairedCor(true.u,true.v,pred.u,pred.v)
  corr.ret <- list("cb_corr_u"=corr.based$corr_u, "cb_corr_v"=corr.based$corr_v,
                   "cb_kendall_u"=corr.based$kendall_u, "cb_kendall_v"=corr.based$kendall_v,
                   "cb_rrmse_u"=rrmse(globalScale(corr.based$scnd_u),globalScale(corr.based$lead_u)),
                   "cb_rrmse_v"=rrmse(globalScale(corr.based$scnd_v),globalScale(corr.based$lead_v)))
  return(list("numeric_dat" = c(procrust.dat$metrics, corr.ret), 
              "pr_matrices" =procrust.dat$paired_matrices, "cb_matrices"=corr.based))

}


globalScale <- function(mat)
{
  test <- matrix(rnorm(100), ncol=10)
  (test - mean(test))/sd(test)
}
matrixDimsAlign <- function(a,b)
{
  if(all(dim(a) != dim(b))) {return(FALSE)}
  TRUE
}

#A u-first orientation.
procrustesPairedUV <- function(fg.u,fg.v,uk.u,uk.v)
{
  if(!(matrixDimsAlign(fg.u, uk.u) & matrixDimsAlign(fg.v, uk.v))) {
    warning("Matrix dimensions aren't aligned")
    return(NA)
  }
  #Enforce order - we want fg first. Given with respect to FG in all cases.
  #My way of doing it
  procrust.paired <- twoMatrixProcrustes(fg.u, uk.u,fg.v,uk.v, scale = FALSE)
  procrust <- vegan::procrustes(fg.u,uk.u,scale=TRUE) #Did we want symmetric or not?
  
  remake <-scale(uk.u,scale=FALSE) %*% procrust$rotation *  procrust$scale
  stopifnot(all(remake == procrust$Yrot) | all(is.nan(remake)))
  if(all(is.nan(procrust$Yrot)))
  {
    message("All NaN case making 0")
    procrust$Yrot[is.nan(procrust$Yrot)] <- 0
  }
  rrmse.u <- rrmse(procrust$Yrot, procrust$X) #ORDER IS PRED, TRUE
  
  ### Three options to try, since its not clear to me how to make the mean centering thing okay.
  #Center along K, center along m, don't center at all
  scaled.uk.m <- scale(t(uk.v),scale=FALSE) #center along elements of m;
  ref.fg.m <- scale(t(fg.v),scale=FALSE) #center along elements of m; 
  
  scaled.uk.k <- t(scale(uk.v,scale=FALSE)) #center along elements of k, this is what I suspect is correct
  ref.fg.k <- t(scale(fg.v,scale=FALSE)) #center along elements of m; 
  if(is.nan(procrust$scale) |  procrust$scale == 0)
  {
	 message("Empty matrix")
    #all 0 matrix, not really sure what to do with this...
    #If NA'd, its a 0 matrix.
    transformed.uk.m <- scaled.uk.m * 0
    transformed.uk.none <- scaled.uk.m * 0
    transformed.uk.k <- scaled.uk.m * 0
  }else
  {
    transformed.uk.k <- solve(procrust$rotation *  procrust$scale) %*% scaled.uk.m
    transformed.uk.m <- solve(procrust$rotation *  procrust$scale) %*% scaled.uk.k
    transformed.uk.none <- solve(procrust$rotation *  procrust$scale) %*% t(uk.v)
  }
  
  rrmse.v.m <- rrmse(transformed.uk.m,ref.fg.m) #centering along m seems to lead out, which doesn't make sense
  rrmse.v.k <- rrmse(transformed.uk.k,ref.fg.k)
  #rrmse.v.none <- rrmse(transformed.uk.none,t(fg.v))
  #stopifnot(rrmse.v.m < rrmse.v.k)
  norm.v.m <- norm(transformed.uk.m-ref.fg.m,"F") #centering along m seems to lead out, which doesn't make sense
  norm.v.k <- norm(transformed.uk.k-ref.fg.k,"F") 
  #norm.v.none <- norm(transformed.uk.none-t(fg.v),"F") 
  #stopifnot(norm.v.m < norm.v.k)
  rrmse.v <- rrmse.v.m
  #correlation calcs
  u.calcs <- calcCorrelationMetrics(procrust$Yrot,procrust$X)
  v.calcs <- calcCorrelationMetrics(ref.fg.m,transformed.uk.m)
  
  #Alternative  procrust.paired version

  rrmse.u.paired <- rrmse(procrust.paired$Yrot_1,procrust.paired$X_1) 
  rrmse.v.paired <- rrmse(procrust.paired$Yrot_2,procrust.paired$X_2)
  u.calcs.paired <- calcCorrelationMetrics(procrust.paired$Yrot_1,procrust.paired$X_1)
  v.calcs.paired <- calcCorrelationMetrics(procrust.paired$Yrot_2,procrust.paired$X_2)

  list("metrics" = list("rmse_v"=rrmse.v, "rmse_u"=rrmse.u, 
       "pearson_v" = v.calcs$pearson, "pearson_u"= u.calcs$pearson, 
       "kendall_v"=v.calcs$kendall, "kendall_u"=u.calcs$kendall,
       "rmse_v_paired"=rrmse.v.paired,"rmse_u_paired"=rrmse.u.paired,
       "pearson_v_paired" = v.calcs.paired$pearson, "pearson_u_paired"= u.calcs.paired$pearson, 
       "kendall_v_paired"=v.calcs.paired$kendall, "kendall_u_paired"=u.calcs.paired$kendall),
       "paired_matrices"=list("lead_v"=procrust.paired$X_2, "lead_u"=procrust.paired$X_1, 
                              "scnd_v"=procrust.paired$Yrot_2, "scnd_u"=procrust.paired$Yrot_1, "swap"=FALSE ))
}

calcCorrelationMetrics <- function(a,b)
{
  list("pearson"=stackAndAssessCorr(a, b),
       "kendall"=stackAndAssessCorr(a, b,cor.type = "kendall"))
}


###########################
#### Greedy correlation calculation and test cases
##########################

###### Test cases ###############
#These were generated by chat GPT on Augus 19, 2024
# Load necessary libraries if not already loaded
if(!requireNamespace("pcaPP", quietly = TRUE)) install.packages("pcaPP")

# Assuming all necessary functions are sourced or available in the environment.

# Helper function to run a test case
run_test <- function(test_name, true.v, fp, verbose = FALSE, cor.type = "pearson") {
  cat("Running Test:", test_name, "\n")
  tryCatch({
    result <- greedyMaxCorr(true.v, fp, verbose = verbose, cor.type = cor.type)
    print(result)
  }, error = function(e) {
    cat("Test failed with error:", e$message, "\n")
  })
  cat("\n")
}
testGreedyMaxCorrEdge <- function()
{
  # Test 1: Empty matrix input
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  run_test("Empty Matrix Input", empty_matrix, empty_matrix)
  
  # Test 2: Matrix with one column input
  one_column_matrix_true <- matrix(1:5, nrow = 5, ncol = 1)
  one_column_matrix_fp <- matrix(5:1, nrow = 5, ncol = 1)
  run_test("One Column Matrix Input", one_column_matrix_true, one_column_matrix_fp)
  
  # Test 3: All-zero correlation matrix input
  # Here true.v and fp are chosen such that all their columns are orthogonal to each other, leading to zero correlation
  all_zero_corr_true <- matrix(c(1,0,0, 0,1,0, 0,0,1), nrow = 3, byrow = TRUE)
  all_zero_corr_fp <- matrix(c(0,1,0, 0,0,1, 1,0,0), nrow = 3, byrow = TRUE)
  run_test("All-Zero Correlation Matrix Input", all_zero_corr_true, all_zero_corr_fp)
  
  # Test 4: Multiple maximum correlations
  # Here, the first and second columns of true.v match equally well with the first column of fp
  multi_max_corr_true <- matrix(c(1, 1, 0, 0, 1, 1), nrow = 2, byrow = TRUE)
  multi_max_corr_fp <- matrix(c(1, 1), nrow = 2)
  run_test("Multiple Maximum Correlations", multi_max_corr_true, multi_max_corr_fp)
  
  # Test 5: Small size but valid inputs with known outcome
  small_true_v <- matrix(c(2, 0, 3, 4, 1, 6), nrow = 2, byrow = TRUE)
  small_fp <- matrix(c(4, 1, 6, 2, 0, 3), nrow = 2, byrow = TRUE)
  run_test("Small Size but Valid Inputs", small_true_v, small_fp)
  
}
testGreedyMaxCorr <- function()
{
  #Manual example for checking performacne
  A <- matrix(rnorm(50),nrow=10)
  B <- matrix(0,nrow=10, ncol = 5)
  for(i in 1:5)
  {
    print(i%%5+1)
    B[,(i%%5+1)] <- A[,i]*(i/3) + rnorm(10,sd=0.1) + rnorm(10,sd=0.03)
    
  }
  cor(A,B)
  #Correct order: 
  test <- greedyMaxCorr(A, B, verbose = TRUE, cor.type = "pearson")
  test.tab <- data.frame("A_order"= test$order.true, "B_order" =test$order.pred) %>% arrange(`A_order`)
  stopifnot(all(test.tab$B_order==c(2,3,4,5,1)))
  
  A.c <- A
  A.c[,1] <- A[,2]
  A.c[,2] <- A[,5]
  A.c[,3] <- A[,4]
  A.c[,4] <- A[,3]
  A.c[,5] <- A[,1]
  test2 <- greedyMaxCorr(A.c, B, verbose = TRUE, cor.type = "pearson")
  test2.tab <- data.frame("A_order"= test2$order.true, "B_order" =test2$order.pred) %>% arrange(`A_order`)
  stopifnot(all(test2.tab$B_order==c(3,1,5,4,2)))
  #Second test- add some zeroes in there.
  B[,3] <- rep(0,10)
  #Correct order: 
  test <- greedyMaxCorr(A, B, verbose = TRUE, cor.type = "pearson")
  test.tab <- data.frame("A_order"= test$order.true, "B_order" =test$order.pred) %>% arrange(`A_order`)
  stopifnot(all(test.tab$B_order==c(2,3,4,5,1)))
  
  #Shuffle A a bit up
  
  
}
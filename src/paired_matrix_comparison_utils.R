GLOBAL_THRESH=7
source("/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/src/evaluateSimR2.R")
`%>%` <- magrittr::`%>%`

#'
#' @param A 
#' @param B 
#'
#' @returns list with objects
#' "lead"- the larger matrix
#' "second"- the smaller matrix, filled with 0s
#' "swap"- if we had to reorder from input
#' NOTE- swapping feature has been removed- always set to FALSE.
#' @export
#'
#' @examples
tidyMatrices <- function(A,B)
{
  SWAP = FALSE #refers to sizing order to get to the same size.
  lead <- as.matrix(A); scnd <- as.matrix(B)
  if(ncol(B) > ncol(A))
  {
    #SWAP = TRUE
    ret <- fillWithZeros(B, A)
    lead <- as.matrix(ret$fp)
  }
  if(nrow(A) != nrow(B))
  {
    warning("Matrices don't have the same number of entries, beware.")
    lead <- emptyCheck(lead, scnd)
    scnd <- emptyCheck(scnd, lead)
  }
  if(ncol(A) > ncol(B))
  {
    #SWAP = TRUE
    ret <- fillWithZeros(A, B)
    scnd <- ret$fp
  }
  
  list("lead"=lead, "second"=scnd, "swap"=SWAP)
}

#Helper function in case either matrix is empty.
emptyCheck <- function(repair.me, reference)
{
  if(nrow(repair.me) == 1 & is.na(repair.me[1,1]))
  {
    repair.me <- matrix(0,ncol =ncol(reference), nrow=nrow(reference))
  }
  repair.me
}
###A greedy implementation to evaluate factor construction, one matrix at a time.
#The main recursive call- requires
#Correlation matrix between remaining factors
#R2 matrix of remaining factors
#signs- vector of the signs assigned to each column
#pairs-a list where each entry is a tuple of (reference vect, pred vect) that match to each other
#returns signs and pairs at the end.
#recurseGreedyMap(r.mat,r2, c(dir), list(dt))
#    matched.cols <- recurseGreedyMapPaired(pcd$r_mat_u, pcd$r_mat_v,pcd$r2_global, c(pcd$dir), list(pcd$order), TRUE, init.dat)

recurseGreedyMapPaired <- function(cor.u, cor.v, r2.global, signs, pairs, global_at_thresh,base.data, global.thresh=GLOBAL_THRESH)
{
  #Check- have we reduced the number enough to avoid needing to be greedy anymore, but can be optimal?
  if(global_at_thresh)
  {
    remaining.cols <- ncol(cor.u) - length(pairs)
    if(sum( remaining.cols < global.thresh))
    {
      #message("Small enough matrix for global approach, proceeding with that")
      used.cols <- do.call("rbind", pairs) #
      remaining.ref.cols <- (1:nrow(cor.v))[!(1:nrow(cor.v) %in% used.cols[,1])] #rows are reference data
      remaining.pred.cols <- (1:ncol(cor.v))[!(1:ncol(cor.v) %in% used.cols[,2])]#cols are pred data
      last.cols <- pairedMatrixAlignmentByCorrelation(base.data$true.u[,remaining.ref.cols], 
                                  base.data$true.v[,remaining.ref.cols], 
                                  base.data$pred.u[,remaining.pred.cols], 
                                  base.data$pred.v[,remaining.pred.cols])
      new.order <- cbind(remaining.ref.cols,remaining.pred.cols[last.cols$new_order]); colnames(new.order) <- c("row", "col")
      new.pairs <- lapply(1:nrow(new.order), function(i) new.order[i,])
      #Note that signs outputted from this have no bearing on previous signs, but we do need them to align. So could be either *1 signs or *-1 signs
      final.signs <- alignGlobalSigns(last.cols$signs,new.pairs, base.data, signs, pairs)
      #Update the signs and pairs, make sure this is right >_<
      return(list("signs"=c(signs, final.signs), 
                  "pairs" = c(pairs,new.pairs)))
    }

  }
    if(all(r2.global == 0) | ncol(r2.global) == 1)
  {
    return(list("signs"=signs, "pairs" = pairs))
  }

  dt <- which(r2.global == max(r2.global), arr.ind = TRUE)[1,] #just go with the first one...
  dir.u <- sign(cor.u[dt[1], dt[2]]);   dir.v <- sign(cor.v[dt[1], dt[2]]);
  if(dir.u != dir.v)
  {
    message("Correlation entry at level ", length(signs),  " has different signs in U and V. ")
    message("This can occur, for instance, if correlation is ~0 in one matrix and high in the other, such that the lowly-correlated factor has an opposite sign by random chance.")
    message("We enforce the sign of the more strongly correlated factor here")
    message("----------------------------------")
    dir.v <- sign(cor.u[dt[1], dt[2]]); #set as the one in U by default
    if(cor.v[dt[1], dt[2]]^2 > cor.u[dt[1], dt[2]]^2)
    {
      dir.v <- sign(cor.v[dt[1], dt[2]]); dir.u <- sign(cor.v[dt[1], dt[2]])
    }
  }

  pairs[[length(pairs) + 1]] <- dt #this tells us which x element goes with which y element 
  
  cor.u[dt[1],] <- 0; cor.u[,dt[2]] <- 0
  cor.v[dt[1],] <- 0; cor.v[,dt[2]] <- 0
  r2.global[dt[1],] <- 0; r2.global[,dt[2]] <- 0
  #recurseGreedyMap(cor,r2, c(signs, dir.u), pairs)
  recurseGreedyMapPaired(cor.u, cor.v, r2.global, c(signs, dir.u), pairs, global_at_thresh,base.data)

}

alignGlobalSigns <- function(new_signs, new_pairs, init.dat, curr_signs, curr_pairs)
{
  order.list <- 1:ncol(init.dat$true.v)
  all.pairs <- c(curr_pairs,new_pairs)
  true.order <- sapply(all.pairs, function(x) x[1]) 
  pred.order <- sapply(all.pairs, function(x) x[2])
  all.signs.pos <- c(curr_signs, new_signs)
  all.signs.neg <- c(curr_signs, -1 * new_signs)
  lead.v <- as.matrix(init.dat$true.v[,true.order]); scnd.v <- matrixSignsProduct(init.dat$pred.v[,pred.order], all.signs.pos)
  lead.u <- as.matrix(init.dat$true.u[,true.order]); scnd.u <- matrixSignsProduct(init.dat$pred.u[,pred.order], all.signs.pos)
  v.pos.cor <- stackAndAssessCorr(lead.v, scnd.v); u.pos.cor <- stackAndAssessCorr(lead.u, scnd.u)
  #Versus
  scnd.v.n <- matrixSignsProduct(init.dat$pred.v[,pred.order], all.signs.neg)
  scnd.u.n <- matrixSignsProduct(init.dat$pred.u[,pred.order], all.signs.neg)
  v.neg.cor <- stackAndAssessCorr(lead.v, scnd.v.n); u.neg.cor <- stackAndAssessCorr(lead.u, scnd.u.n)

  if((v.neg.cor+u.neg.cor)^2 >(v.pos.cor + u.pos.cor)^2)
  {
    return(-1 * new_signs)
  }
  return(new_signs)
  
  
}


if(FALSE)
{
  true.v <- model.results.sub$finngen$`BIC-dev_K-GRID`$V
  true.u <- model.results.sub$finngen$`BIC-dev_K-GRID`$U
  pred.v <- model.results.sub$ukbb$`BIC-dev_K-GRID`$V
  pred.u <- model.results.sub$ukbb$`BIC-dev_K-GRID`$U
}

#' Reorder pairs matrices to maximize the global correlation with respect to a 
#' true U and V and predicted U and V. Must all have same dimensions, recommend padding with 0s
#'
#' @param true.u 
#' @param pred.u - must be same dimensions as true.u
#' @param true.v 
#' @param pred.v - must be same dimensions as true.v
#' @param cor.type pearson or kendall
#' @param global_at_thresh do you want to use the global optima finding approach when there are sufficiently few columns?
#'        recommended setting is TRUE
#'
#' @return a list containing the global correlation score, the order of the reference columns, the order of the predicted columns, and the accompanying sign changes.
#' @export
#'
#' @examples
greedyMaxPairedCor <- function(true.u, true.v, pred.u, pred.v, cor.type="pearson", global_at_thresh = TRUE)
{
  if(!(matrixDimsAlign(true.u, pred.u) & matrixDimsAlign(true.v, pred.v))) {
    warning("Matrix dimensions aren't aligned. Recommend padding missing cols with 0.")
    return(NA)
  }
  init.dat <- list("true.v"=true.v, "pred.v"=pred.v,
                   "true.u" =true.u, "pred.u" = pred.u)
  lower.limit = calcCorrelationMetrics(true.v, pred.v)
  if(ncol(init.dat$true.v) < GLOBAL_THRESH & global_at_thresh)
  {
    message("Small enough matrix for global approach, proceeding with that")
    short.out <- pairedMatrixAlignmentByCorrelation(init.dat$true.u, init.dat$true.v, init.dat$pred.u, init.dat$pred.v)
    #list("U_r2"=l_r2, "V_r2"=f_r2, "reordered_U"=lp,"reordered_V"= fp, "new_order"=ordering[[order_i]], "signs"=unlist(sign.grid[sign_i,])))
    new.order <- cbind(1:ncol(init.dat$true.v),short.out$new_order); colnames(new.order) <- c("row", "col")
    matched.cols <- list("signs"=short.out$signs, 
                         "pairs" = lapply(1:nrow(new.order), function(i) new.order[i,]))
  }else
  {
    pcd <- prepForRecursiveGreedySearch(init.dat, cor.type)
    matched.cols <- recurseGreedyMapPaired(pcd$r_mat_u, pcd$r_mat_v,pcd$r2_global, c(pcd$dir), list(pcd$order), TRUE, init.dat)
    # function(cor.u, cor.v, r2.global, signs, pairs, global_at_thresh,base.data, global.thresh=GLOBAL_THRESH)
  }

  order.list <- 1:ncol(true.v)
  true.order <- sapply(matched.cols$pairs, function(x) x[1]) %>%  c(., order.list[!(1:ncol(true.v) %in% .)]) #add on any missing columns, should be none
  pred.order <- sapply(matched.cols$pairs, function(x) x[2])%>%  c(., order.list[!(1:ncol(true.v) %in% .)])  #add on any missing columns, should be none
  pred.signs <-  c(matched.cols$signs, rep(1,(ncol(true.v) - length(matched.cols$signs)))) #This adds in 1s for any columns that were unassigned

  lead.v <- as.matrix(true.v[,true.order]); scnd.v <- matrixSignsProduct(pred.v[,pred.order], pred.signs)
  lead.u <- as.matrix(true.u[,true.order]); scnd.u <- matrixSignsProduct(pred.u[,pred.order], pred.signs)
  
  c_v <- stackAndAssessCorr(lead.v, scnd.v)
  c_u <- stackAndAssessCorr(lead.u, scnd.u)
  
    #set the true order to what it was, and the pred order as modified:
    order.f <- order(true.order)
    true.order = true.order[order.f]
    pred.order = pred.order[order.f]
    pred.signs = pred.signs[order.f]
    
    #Some sanitcy checks while under testing....
    if(!all.equal(c_v,stackAndAssessCorr(as.matrix(true.v[,true.order]), matrixSignsProduct(pred.v[,pred.order], pred.signs))))
    {
      message("Look here")
    }
    if(!all.equal(c_u,stackAndAssessCorr(as.matrix(true.u[,true.order]), matrixSignsProduct(pred.u[,pred.order], pred.signs))))
    {
      message("Look here")
    }
    stopifnot(all.equal(c_v,stackAndAssessCorr(as.matrix(true.v[,true.order]), matrixSignsProduct(pred.v[,pred.order], pred.signs))))
    stopifnot(all.equal(c_u,stackAndAssessCorr(as.matrix(true.u[,true.order]), matrixSignsProduct(pred.u[,pred.order], pred.signs))))
    return(list("corr_v"=c_v, "corr_u"=c_u,
                "kendall_v"=stackAndAssessCorr(lead.v, scnd.v,cor.type = "kendall"),"kendall_u"=stackAndAssessCorr(lead.u, scnd.u,cor.type = "kendall"),
                "order.true"=true.order, "order.pred"=pred.order, "signs"=pred.signs,
                "lead_v"=lead.v, "lead_u"=lead.u, "scnd_v"=scnd.v, "scnd_u"=scnd.u, "swap"=FALSE))
    #corr.based$kendall_v
  }
  
  


#' Helper function to prepare data for recursive greedy search
#'
#' @param init.dat contains all of the matrices (ref v, pred v, ref u, pred u)
#' @param cor.type 
#'
#' @returnlist containing the correlation matrices (r_mat_u/v), the global r2 matrix, 
#' the sign of the first choice and the coordinates of the first choice
#' @export
#'
#' @examples
prepForRecursiveGreedySearch <- function(init.dat, cor.type)
{
  #true is on x axis, pred is on y axis
  r.mat.u <- customCorr(init.dat$true.u, init.dat$pred.u, cor.type=cor.type);  r.mat.u[is.na(r.mat.u)] <- 0
  r.mat.v <- customCorr(init.dat$true.v, init.dat$pred.v, cor.type=cor.type);  r.mat.v[is.na(r.mat.v)] <- 0
  r2.u <-  r.mat.u^2;  r2.v <-  r.mat.v^2;  
  #r2.global <- r2.u + r2.v
  #MAJOR CHANGE- doing based on correlation to require directions to be consistent:
  r2.global <- (r.mat.u + r.mat.v)^2
  #Get the top entry
  dt <- which(r2.global == max(r2.global), arr.ind = TRUE)[1,]
  stopifnot(max(r2.global) == r2.global[matrix(dt,ncol=2)])
  #Get the sign of the top entry
  dir.u <- sign(r.mat.u[dt[1], dt[2]]);   dir.v <- sign(r.mat.v[dt[1], dt[2]]);
  if(dir.u != dir.v)
  {
    message("Top correlation entry has different signs in U and V. Review this case manually")
    return(NULL)
  }
  
  #Zero-out the rows and columns corresponding to entries we have already selected.
  r.mat.v[dt[1],] <- 0; r.mat.v[,dt[2]] <- 0
  r.mat.u[dt[1],] <- 0; r.mat.u[,dt[2]] <- 0
  r2.global[dt[1],] <- 0; r2.global[,dt[2]] <- 0
  return(list("r_mat_u" = r.mat.u, "r_mat_v"=r.mat.v, "r2_global"=r2.global, "dir"=dir.u, "order"=dt))
}


#Helper function for pairedMatrixAlignmentByCorrelation- see documentation below
evaluate_error <- function(trueL, trueF, lp, fp)
{
  pairedMatrixAlignmentByCorrelation(trueL, trueF, lp, fp)
}
#' Comphrensive and simple R2 between matrices
#' Based on code from Yuan He, this samples every possible ordering of columns and sign assignment to find the very best one
#' It then returns the R2 for this optimal arrangement.
#' @param trueL - known U matrix
#' @param trueF - known V matrix
#' @param lp - predicted U matrix
#' @param fp - predicted V matrix
#'
#' @return
#' @export
#'
#' @examples
pairedMatrixAlignmentByCorrelation <- function(trueL, trueF, lp, fp){
  #Special case- all entries are 0
  if(all(fp == 0) & all(lp == 0))
  {
    #check that all trueL and trueF aren't zero
    
    if(all(trueL == 0) & all(trueF == 0))
    {
      message("Input is all empty matrices.")
     
    }
    return(list("U_r2"=0, "V_r2"=0, "reordered_U"=lp,"reordered_V"= fp, "new_order"=1:ncol(trueL), "signs"=rep(1,ncol(trueL))))
  }
  if(ncol(fp) < ncol(trueF)){
    dif = ncol(trueF) - ncol(fp)
    fp = cbind(fp, matrix(rep(0, nrow(fp) * dif), ncol = dif))
    lp = cbind(lp, matrix(rep(0, nrow(lp) * dif), ncol = dif))
  }
  rankK = ncol(trueF)
  if(rankK > GLOBAL_THRESH)
  {
    message("There are too many latent factors for comprehensive examination of options. Do a greedy search.")
    return()
  }
  suppressWarnings(library('combinat'))
  ordering = permn(seq(1,rankK))
  #Signs
  #n = 5
  n=rankK
  sign.grid <- expand.grid(replicate(n, c(1,-1), simplify = FALSE))
  #all.combs <- apply(sign.grid,1,function(x) lapply(ordering, function(y) x*y))
  #long.list <- lapply(all.combs, unlist, use.names=FALSE)
  ntest <- length(ordering) * nrow(sign.grid)
  test_i = 1
  f_cor = rep(0, ntest); l_cor = rep(0, length(ntest))
  for(ord in 1:length(ordering)){
    for(i in 1:nrow(sign.grid))
    {
      #matrixSignsProduct(scnd[,procestes.corr$mapping_reorder],procestes.corr$mapping_signs))
      #f_cor[test_i] = customCorr(as.vector(trueF), as.vector(matrixSignsProduct(fp[,ordering[[ord]]],sign.grid[i,])))^2
      #l_cor[test_i] = customCorr(as.vector(trueL), as.vector(matrixSignsProduct(lp[,ordering[[ord]]],sign.grid[i,])))^2
      f_cor[test_i] = customCorr(as.vector(trueF), as.vector(matrixSignsProduct(fp[,ordering[[ord]]],sign.grid[i,])))
      l_cor[test_i] = customCorr(as.vector(trueL), as.vector(matrixSignsProduct(lp[,ordering[[ord]]],sign.grid[i,])))
      test_i = test_i+1
    }
  }
  
  #All sign options:
  ord_sum = (f_cor + l_cor)^2
  opt_i = which.max(ord_sum)
  if(length(opt_i) > 1)
  {
    message("Multiple optimal orientations, selecting the first....")
    opt_i = opt_i[1]
  }
  #Translate the global index into the order/sign combination.
  order_i <- ceiling(opt_i / nrow(sign.grid))
  sign_i <- opt_i %% nrow(sign.grid) #Only exception is if 
  #Test
  #best.cor <- customCorr(as.vector(trueF), as.vectfor(fp[,ordering[[order_i]]] %*% diag(sign.grid[sign_i,])))^2 + customCorr(as.vector(trueL), as.vector(lp[,ordering[[order_i]]] %*% diag(sign.grid[sign_i,])))^2
  best.cor <- (customCorr(as.vector(trueF), as.vector(fp[,ordering[[order_i]]] %*% diag(sign.grid[sign_i,]))) + customCorr(as.vector(trueL), as.vector(lp[,ordering[[order_i]]] %*% diag(sign.grid[sign_i,]))))^2
  stopifnot(best.cor == max(ord_sum) )
  lp = lp[,ordering[[order_i]]] %*% diag(sign.grid[sign_i,])
  fp = fp[,ordering[[order_i]]] %*% diag(sign.grid[sign_i,])
  
  l_r2 = customCorr(as.vector(trueL), as.vector(lp))^2
  f_r2 = customCorr(as.vector(trueF), as.vector(fp))^2
  
  #new_order is the order in which the predicted matrix should be to maximize global correlation
  #signs is the sign of which those should be to maximize global correlation
  return(list("U_r2"=l_r2, "V_r2"=f_r2, "reordered_U"=lp,"reordered_V"= fp, "new_order"=ordering[[order_i]], "signs"=unlist(sign.grid[sign_i,])))
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




#############Test cases:
### Written with the help of ChatGPT, 8/20/2024
# Unit test for empty matrices
test_empty_matrices <- function() {
  true.u <- matrix(numeric(0), nrow = 0, ncol = 0)
  pred.u <- matrix(numeric(0), nrow = 0, ncol = 0)
  true.v <- matrix(numeric(0), nrow = 0, ncol = 0)
  pred.v <- matrix(numeric(0), nrow = 0, ncol = 0)
  
  result <- pairedMatrixAlignmentByCorrelation(true.u, true.v, pred.u, pred.v)
  procrust.result <- twoMatrixProcrustes(true.u, pred.u,true.v,pred.v)
  
  stopifnot(all(result$U_r2 == 0))
  stopifnot(all(result$V_r2 == 0))
  stopifnot(all(dim(result$reordered_U) == c(0, 0)))
  stopifnot(all(dim(result$reordered_V) == c(0, 0)))
  
  stopifnot(all(procrust.result$Yrot_1 == 0))
  stopifnot(all(procrust.result$Yrot_2 == 0))
  stopifnot(procrust.result$rotation == 1)

  message("Empty matrices test passed!")
}

# Unit test for matrices with one column
test_one_column_matrices <- function() {
  true.u <- matrix(rnorm(20), nrow = 20, ncol = 1)
  pred.u <- true.u + rnorm(20, sd = 0.5)
  true.v <- matrix(rnorm(10), nrow = 10, ncol = 1)
  pred.v <- true.v + rnorm(10, sd = 0.5)
  
  result <- evaluate_error(true.u, true.v, pred.u, pred.v)
  procrust.result <- twoMatrixProcrustes(-1*true.u, pred.u,-1*true.v,pred.v)
  procrust.result.pos <- twoMatrixProcrustes(true.u, pred.u,true.v,pred.v)
  stopifnot(procrust.result$rotation == -1)
  stopifnot(procrust.result.pos$rotation == 1)
  stopifnot(result$U_r2 > 0)
  stopifnot(result$V_r2 > 0)
  stopifnot(all(dim(result$reordered_U) == dim(pred.u)))
  stopifnot(all(dim(result$reordered_V) == dim(pred.v)))
  message("One column matrices test passed!")
}

# Unit test for all-zero correlation matrices
test_all_zero_correlation <- function() {
  true.u <- matrix(0, nrow = 20, ncol = 5)
  pred.u <- matrix(0, nrow = 20, ncol = 5)
  true.v <- matrix(0, nrow = 10, ncol = 5)
  pred.v <- matrix(0, nrow = 10, ncol = 5)
  
  result <- evaluate_error(true.u, true.v, pred.u, pred.v)
  procrust.result <- twoMatrixProcrustes(true.u, pred.u,true.v,pred.v)
  stopifnot(result$U_r2 == 0)
  stopifnot(result$V_r2 == 0)
  stopifnot(all(dim(result$reordered_U) == dim(pred.u)))
  stopifnot(all(dim(result$reordered_V) == dim(pred.v)))
  stopifnot(procrust.result$rotation == 1)
  message("All-zero correlation matrices test passed!")
}

########################################################################
### Correlation-based method specific tests
########################################################################
testEvaluateError <- function()
{
  set.seed(2)
  #Simulate U
  true.u <- matrix(rnorm(100), nrow=20, ncol=5)
  mod.u <- true.u + matrix(rnorm(100,sd = 0.5), nrow=20, ncol=5)
  ref.r2.u <- cor(as.vector(true.u), as.vector(mod.u))^2
  #Simulate V
  true.v <- matrix(rnorm(50), nrow=10, ncol=5)
  mod.v <- true.v + matrix(rnorm(50,sd = 0.5), nrow=10, ncol=5)
  ref.r2.v <- cor(as.vector(true.v), as.vector(mod.v))^2
  
  #move things around
  cp.u <- mod.u
  cp.u[,1] <- mod.u[,2] * -1
  cp.u[,2] <- mod.u[,1] * -1
  
  cp.v <- mod.v
  cp.v[,1] <- mod.v[,2] * -1
  cp.v[,2] <- mod.v[,1] * -1
  
  check.dat <- evaluate_error(true.u, true.v, cp.u, cp.v)
  stopifnot(check.dat$U_r2 == ref.r2.u)
  stopifnot(check.dat$V_r2 == ref.r2.v)
  stopifnot(check.dat$reordered_U == mod.u)
  stopifnot(check.dat$reordered_V == mod.v)
  message("Sample test passed!")
}

#edge case I am concerned with- what if U and v map to different places? Does it handle that okay?
testEvaluateGreedy <- function()
{
  set.seed(3)
  #Simulate U
  #True- 11-11
  true.u <- cbind(matrix(rnorm(1000), nrow=100, ncol=10), c(1,rep(0,99)))
  mod.u <- true.u + matrix(rnorm(1100,sd = 0.5), nrow=100, ncol=11) #added noise
  ref.r2.u <- cor(as.vector(true.u), as.vector(mod.u))^2
  #Simulate V
  true.v <- matrix(rnorm(200), nrow=20, ncol=10)
  true.v <- cbind(true.v,true.v[,10]*0.4+rnorm(20,sd=0.1)) #v's column 11 could map to 10
  
  mod.v <- true.v + matrix(rnorm(220,sd = 0.5), nrow=20, ncol=11)
  ref.r2.v <- cor(as.vector(true.v), as.vector(mod.v))^2
  
  #move things around
  cp.u <- mod.u
  cp.u[,1] <- mod.u[,2] * -1
  cp.u[,2] <- mod.u[,1] * -1
  cp.u[,11] <- mod.u[,9] 
  cp.u[,9] <- mod.u[,11]
  
  cp.v <- mod.v
  cp.v[,1] <- mod.v[,2] * -1
  cp.v[,2] <- mod.v[,1] * -1
  cp.v[,11] <- mod.v[,9] 
  cp.v[,9] <- mod.v[,11]
  
  #true.u, pred.u, true.v, pred.v
  check.dat <- greedyMaxPairedCor(true.u, true.v, cp.u, cp.v)
  stopifnot(check.dat$corr_u^2 == ref.r2.u)
  stopifnot(check.dat$corr_v^2 == ref.r2.v)
  stopifnot(all(matrixSignsProduct(cp.u[,check.dat$order.pred],check.dat$signs)  == mod.u))
  stopifnot(all(matrixSignsProduct(cp.v[,check.dat$order.pred],check.dat$signs)  == mod.v))
  message("More complex test passed!")
}

########################################################################
### Paired procrustes-based method specific tests
#######################################################################

testProcrustIdentical <- function()
{
    set.seed(3)
    #Simulate U
    true.u <- matrix(rnorm(100), nrow=20, ncol=5)
    mod.u <- true.u + matrix(rnorm(100,sd = 0.01), nrow=20, ncol=5)
    
    #Simulate V to be the exact same
    true.v <- true.u
    mod.v <- mod.u
  
    procrust.id <- twoMatrixProcrustes(true.u, mod.u,true.v,mod.v)
    
    #expected outcome: no rotation, they are in the right order to begin with..
    stopifnot(all(diag(round(procrust.id$rotation, digits = 2)) == 1))
    message("Identical matrices procrustes test passed!")
  }

testProcrustIdenticalSwappedCol <- function()
{
  set.seed(3)
  #Simulate U
  true.u <- matrix(rnorm(100), nrow=20, ncol=5)
  mod.u <- true.u + matrix(rnorm(100,sd = 0.01), nrow=20, ncol=5)
  u.pred <- mod.u
  u.pred[,2] <- mod.u[,1]
  u.pred[,1] <- -1 * mod.u[,2]
  #Simulate V to be the exact same
  true.v <- true.u
  mod.v <- mod.u
  v.pred <- mod.v
  v.pred[,2] <- mod.v[,1]
  v.pred[,1] <- -1 * mod.v[,2]
  procrust.id <- twoMatrixProcrustes(true.u, u.pred,true.v,v.pred)
  stopifnot(round(procrust.id$rotation[2,1],digits = 2) == 1)
  stopifnot(round(procrust.id$rotation[1,2],digits = 2) == -1)
  stopifnot(all(round(diag(procrust.id$rotation),digits = 2) == c(0,0,1,1,1)))
  message("Identical matrices w/ swapped columns procrustes test passed!")
}

#test where teh column orders are correct, but U and V are different
testProcrustMatched <- function()
{
  set.seed(5)
  #Simulate U
  true.u <- matrix(rnorm(100), nrow=20, ncol=5)
  pred.u <- true.u + matrix(rnorm(100,sd = 0.1), nrow=20, ncol=5)
  #Simulate V to be the exact same
  true.v <- matrix(rnorm(35), nrow = 7, ncol = 5)
  pred.v <- true.v + rnorm(35, sd = 0.1)
  procrust.id <- twoMatrixProcrustes(true.u, pred.u,true.v,pred.v)
  stopifnot(all(round(diag(procrust.id$rotation),digits = 2) == 1))
  message("Distinct matrices w/ ordered columns procrustes test passed!")
}

#U and V are different, column orders are different.
testProcrustMatchedSwapped <- function()
{
  set.seed(6)
  #Simulate U
  true.u <- matrix(rnorm(100), nrow=20, ncol=5)
  mod.u <- true.u + matrix(rnorm(100,sd = 0.1), nrow=20, ncol=5)
  u.pred <- mod.u
  u.pred[,2] <- mod.u[,1]
  u.pred[,1] <- -1 * mod.u[,5]
  u.pred[,5] <- mod.u[,2]
  #Simulate V to be the exact same
  true.v <- matrix(rnorm(35), nrow = 7, ncol = 5)
  mod.v <- true.v + rnorm(35, sd = 0.1)
  v.pred <- mod.v
  v.pred[,2] <- mod.v[,1]
  v.pred[,1] <- -1 * mod.v[,5]
  v.pred[,5] <- mod.v[,2]

  #X axis is in pred matrix coordins, y is in true
  procrust.id <- twoMatrixProcrustes(true.u, u.pred,true.v,v.pred)
  stopifnot(round(procrust.id$rotation[2,1],digits = 2) == 1)
  stopifnot(round(procrust.id$rotation[1,5],digits = 2) == -1)
  stopifnot(round(procrust.id$rotation[5,2],digits = 2) == 1)
  stopifnot(all(round(diag(procrust.id$rotation),digits = 1) == c(0,0,1,1,0)))
  message("Distinct matrices w/ disordered columns procrustes test passed!")
}

#Are we truly minmizing the distance?
testProcrustVSVegan <- function()
{
  set.seed(6)
  #Simulate U
  true.u <- matrix(rnorm(100), nrow=20, ncol=5)
  mod.u <- true.u + matrix(rnorm(100,sd = 0.1), nrow=20, ncol=5)
  u.pred <- mod.u
  u.pred[,2] <- mod.u[,1]
  u.pred[,1] <- -1 * mod.u[,5]
  u.pred[,5] <- mod.u[,2]
  #Simulate V to be the exact same
  true.v <- matrix(rnorm(35), nrow = 7, ncol = 5)
  mod.v <- true.v + rnorm(35, sd = 0.1)
  v.pred <- mod.v
  v.pred[,2] <- mod.v[,1]
  v.pred[,1] <- -1 * mod.v[,5]
  v.pred[,5] <- mod.v[,2]
  
  #X axis is in pred matrix coordins, y is in true
  procrust.id <- twoMatrixProcrustes(true.u, u.pred,true.v,v.pred)
  vegan.style <- vegan::procrustes(true.u,u.pred, scale = FALSE)

  procrust.rotate <- t(solve(vegan.style$rotation *  vegan.style$scale) %*% t(v.pred))
  vegan.sum <- rrmse(scale(vegan.style$Yrot), scale(true.u)) + rrmse(scale(procrust.rotate), scale(true.v))
  
  paired.sum <-   rrmse(scale(procrust.id$Yrot_1), scale(true.u)) +   rrmse(scale(procrust.id$Yrot_2), scale(true.v))
  
  stopifnot(vegan.sum > paired.sum)

  stopifnot(all(round(vegan.style$rotation, digits=1) ==   round(procrust.id$rotation, digits=1)))
  message("Paired procrustes matches published procrustes in easy case.")
}
testProcrustVSVeganHard <- function()
{
  set.seed(7)
  #Simulate U, which is noisy
  true.u <- matrix(rnorm(100), nrow=20, ncol=5)
  mod.u <- true.u + matrix(rnorm(100,sd = 1), nrow=20, ncol=5)
  u.pred <- mod.u
  u.pred[,2] <- mod.u[,1]
  u.pred[,1] <- -1 * mod.u[,5]
  u.pred[,5] <- mod.u[,2]
  #Simulate V to be the exact same, but less noisy
  true.v <- matrix(rnorm(35), nrow = 7, ncol = 5)
  mod.v <- true.v + rnorm(35, sd = 0.1)
  v.pred <- mod.v
  v.pred[,2] <- mod.v[,1]
  v.pred[,1] <- -1 * mod.v[,5]
  v.pred[,5] <- mod.v[,2]
  #We expect there to be a more substantial boost
  #X axis is in pred matrix coordins, y is in true
  procrust.id <- twoMatrixProcrustes(true.u, u.pred,true.v,v.pred)
  vegan.style <- vegan::procrustes(true.u,u.pred, scale = FALSE)
  
  procrust.rotate <- t(solve(vegan.style$rotation *  vegan.style$scale) %*% t(v.pred))
  vegan.sum <- rrmse(scale(vegan.style$Yrot), scale(true.u)) + rrmse(scale(procrust.rotate), scale(true.v))
  paired.sum <-   rrmse(scale(procrust.id$Yrot_1), scale(true.u)) +   rrmse(scale(procrust.id$Yrot_2), scale(true.v))
  
  stopifnot(vegan.sum > paired.sum)
  
  stopifnot(all(round(vegan.style$rotation, digits=0) ==   round(procrust.id$rotation, digits=0)))
  message("Paired procrustes matches published procrustes in hard case.")
}

testProcrustVSVeganExperimental <- function()
{
  set.seed(8)
  #Simulate U, which is noisy
  true.u <- matrix(rnorm(100), nrow=20, ncol=5)
  mod.u <- true.u + matrix(rnorm(100,sd = 1), nrow=20, ncol=5)
  u.pred <- mod.u

  #Simulate V to have a single swap in there
  true.v <- matrix(rnorm(35), nrow = 7, ncol = 5)
  mod.v <- true.v + rnorm(35, sd = 0.1)
  v.pred <- mod.v
  v.pred[,2] <- mod.v[,1]
  v.pred[,1] <- mod.v[,2]
  #We expect there to be a more substantial boost
  #X axis is in pred matrix coordins, y is in true
  procrust.id <- twoMatrixProcrustes(true.u, u.pred,true.v,v.pred)
  vegan.style <- vegan::procrustes(true.u,u.pred, scale = FALSE)
  
  procrust.rotate <- t(solve(vegan.style$rotation *  vegan.style$scale) %*% t(v.pred))
  vegan.sum <- rrmse(scale(vegan.style$Yrot), scale(true.u)) + rrmse(scale(procrust.rotate), scale(true.v))
  paired.sum <-   rrmse(scale(procrust.id$Yrot_1), scale(true.u)) +   rrmse(scale(procrust.id$Yrot_2), scale(true.v))
  
  stopifnot(vegan.sum > paired.sum)
    message("Paired procrustes matches published procrustes in mismatched case.")
}

#suggested by chatGPT, and updated to match my aims
testProcrustHighNoise <- function() {
  set.seed(10)
  true.u <- matrix(rnorm(100), nrow=20, ncol=5)
  mod.u <- true.u + matrix(rnorm(100, sd=5), nrow=20, ncol=5)
  true.v <- matrix(rnorm(35), nrow=7, ncol=5)
  mod.v <- true.v + matrix(rnorm(35, sd=5), nrow=7, ncol=5)
  
  procrust.id <- twoMatrixProcrustes(true.u, mod.u, true.v, mod.v)
  stopifnot(all(!is.na(procrust.id$rotation)))
  #Make sure it actually helps
  with.rot  <- rrmse(scale(procrust.id$Yrot_1), scale(procrust.id$X_1)) +  rrmse(scale(procrust.id$Yrot_2), scale(procrust.id$X_2)) 
  no.rot  <- rrmse(scale(true.u), scale(mod.u)) +  rrmse(scale(true.v), scale(mod.v)) 
  stopifnot(with.rot < no.rot)
  message("High noise procrustes test passed!")
}



# Running the tests
testPairedMappingMethods <- function()
{
  test_empty_matrices()
  test_one_column_matrices()
  test_all_zero_correlation()
  testEvaluateError()  # Original function test
  testEvaluateGreedy()
  
  testProcrustMatched()
  testProcrustIdenticalSwappedCol()
  testProcrustIdentical()
  
  testProcrustMatchedSwapped()
  testProcrustVSVeganHard()
  testProcrustVSVeganExperimental()
  testProcrustVSVegan()
  testProcrustHighNoise() #Tests performance when matrices are uncorrelated, and that our rotation is better than no rotation at all.
  message("All tests passed!")
}

####

########


#Procrustes-style approach
#This was based on the code from vegan, updated for my purposes here.

twoMatrixProcrustes <- function (X_1, Y_1,X_2,Y_2, scale = FALSE, symmetric = FALSE,  scores = "sites", ...)
  {
    if (nrow(X_1) != nrow(Y_1))
      stop(gettextf("matrices have different number of rows: %d and %d",
                    nrow(X_1), nrow(Y_1)))
    if (nrow(X_2) != nrow(Y_2))
      stop(gettextf("matrices have different number of rows: %d and %d",
                    nrow(X_2), nrow(Y_2)))
    if (ncol(X_1) < ncol(Y_1)) {
      warning("X has fewer axes than Y: X adjusted to comform Y\n")
      addcols <- ncol(Y_1) - ncol(X_1)
      for (i in 1:addcols) X <- cbind(X_1, 0)
    }
    if(all(Y_1 == 0) & all(Y_2 == 0))
    {
      warning("All entries in target matrices are 0.")
      reslt <- list(Yrot_1 = Y_1, X_1 = (X_1),
                    Yrot_2 = Y_2, X_2 = (X_2), ss = NULL, rotation = 1,
                    translation = NULL, 
                    scale_1 = 1,scale_2=1, xmean_1 = apply(X_1, 2, mean),
                    xmean_2 = apply(X_2, 2, mean),
                    symmetric = symmetric, call = match.call())
      return(reslt)
    }
    ctrace <- function(MAT) sum(MAT^2) #define
    c_2 <- 1; c_1 <- 1
    if (symmetric) {
      stop("not implemented)")
      X <- scale(X, scale = FALSE)
      Y <- scale(Y, scale = FALSE)
      X <- X/sqrt(ctrace(X))
      Y <- Y/sqrt(ctrace(Y))
    }

    xmean.1 <- apply(X_1, 2, mean)
    ymean.1 <- apply(Y_1, 2, mean)
    xmean.2 <- apply(X_2, 2, mean)
    ymean.2 <- apply(Y_2, 2, mean)
    if (!symmetric) {
      #Based on https://rpubs.com/mengxu/procrustes_analysis, and in order to maintain consistency, need to scale entire matrices:
      #Not sure if this makes it symmetric or not though
      X_1 <- preProcrustScale(X_1)
      Y_1 <- preProcrustScale(Y_1)
      X_2 <- preProcrustScale(X_2)
      Y_2 <- preProcrustScale(Y_2)
      #X_1 <- scale(X_1, scale = FALSE)
      #Y_1 <- scale(Y_1, scale = FALSE)
      #X_2 <- scale(X_2, scale = FALSE)
      #Y_2 <- scale(Y_2, scale = FALSE)

    }
    XY_1 <- crossprod(X_1, Y_1)
    XY_2 <- crossprod(X_2, Y_2)
    sol <- svd(XY_1+XY_2)
    A <- sol$v %*% t(sol$u)
    if (scale) {
      stop("not implemented")
      c_1 <- sum(sol$d)/ctrace(Y_1)
      c_2 <- sum(sol$d)/ctrace(Y_2)
    }
    Yrot_1 <- c_1 * Y_1 %*% A
    Yrot_2 <- c_2 * Y_2 %*% A
    ## Translation (b) needs scale (c) although Mardia et al. do not
    ## have this. Reported by Christian Dudel.
    #b <- xmean - c * ymean %*% A
    #R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)
    reslt <- list(Yrot_1 = Yrot_1, X_1 = X_1,
                  Yrot_2 = Yrot_2, X_2 = X_2, ss = NULL, rotation = A,
                  translation = NULL, 
                  scale_1 = c_1,scale_2=c_2, xmean_1 = xmean.1,
                  xmean_2 = xmean.2,
                  symmetric = symmetric, call = match.call())
    reslt$svd <- sol
    class(reslt) <- "procrustes"
    reslt
  }


preProcrustScale <- function(X)
{
  X_ret <- scale(X, scale = FALSE)
  sn <- norm(X_ret, type = "F")
  if(sn == 0 &  all(X==0))
  {
    matrix(0, nrow=nrow(X), ncol =ncol(X))
  }else
  {
    X_ret/(sn / (ncol(X_ret) * nrow(X_ret)))
  } 
}

#A few test cases:
#what happens when the matrix is all 0s
#What happens when the matrix is empty
#Two identical matrices
#
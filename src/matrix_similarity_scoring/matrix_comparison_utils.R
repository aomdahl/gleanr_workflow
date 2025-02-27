#Condensed code extending functionality explored in `comparison_methods.Rmd`
#The purpose of this script is to compare 2 GWAS factorizations for similarity, using standard
#ML-style accuracy scores.


#setup
#source("~/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc//src/plot_functions.R")
source("~/scratch16-abattle4/ashton/snp_networks/custom_l1_factorization/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/matrix_similarity_scoring/evaluateSimR2.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/matrix_similarity_scoring/paired_matrix_comparison_utils.R")
pacman::p_load(ggplot2, magrittr, tidyr)

########Helper functions


#'Get the correlation via procustes, standard CARET, PR metrics of 2 matrices
#'
#' @param A the "reference matrix"
#' @param B the
#'
#' @return list object containing matrices with their columns ordered, mapped via procstes analysis, CARET package reports, PR curve metrics
#' @export

#For testing:
#A <- matrix(rnorm(1000), nrow = 100)
#B <- matrix(rnorm(1000), nrow = 100)

alignAndCompare <- function(ret_obj1, ret_obj2)
{
  v.dat <- prepMatricesForAnalysis(ret_obj1$V,ret_obj2$V); 
  u.dat <- prepMatricesForAnalysis(ret_obj1$U,ret_obj2$U)
  stopifnot(!v.dat$swap)
  fg.v <- v.dat$lead; uk.v <- v.dat$second
  fg.u <- u.dat$lead; uk.u <- u.dat$second
  #Get the matchup
  svd.version <- greedyMaxPairedCor(fg.u,fg.v, uk.u,uk.v)
  #align UKBB to match with Finngen
  uk.u.new = matrixSignsProduct(uk.u[,svd.version$order.pred], svd.version$signs)
  uk.v.new = matrixSignsProduct(uk.v[,svd.version$order.pred], svd.version$signs)
  list("adj.v"=cor.test(y=c(uk.v.new), x=c(fg.v)),"adj.u"= cor.test(y=c(uk.u.new), x=c(fg.u)))
}

#' Align matrices A and B via correlation and via procrustes
#'
#' @param A simple matrix
#' @param B simple matrix
#' @param corr.type which type to guide scoring by.
#'
#' @return an object containing
#'  cb_results: a bunch of data from the correlation-based analysis, including the order
#'  cb_matrices: the output matrices based on this ordering
#'  pr_results: a bunch of data from the procrustes-based
#'  pr_matrices:the output matrices based on this ordering (the raw matrices, reordered and signed, no other scaling)
#' @export
#'
#' @examples
alignMatrices <- function(A,B, corr.type = "pearson")
{
  #First, align the matrices to have the same dimensions
  library(PRROC)
  prepped <- prepMatricesForAnalysis(A,B)
  lead <- prepped$lead;  scnd <- prepped$second; SWAP <- prepped$swap

  #Use my homeade correlation matching
  pseudo.procrust.results <- corrBasedScoring(lead, scnd, corr.method = corr.type)
  pseudo.ordered.matrices <- pseudo.procrust.results$corr_ordered_matrices

  #Use a published procrustes package too
  procrust.results <- procrustesBasedScoring(lead,scnd, SWAP)
  procrust.ordered.matrices <- procrust.results$procrust_ordered_matrices

  list("cb_results"=pseudo.procrust.results, "cb_matrices"=pseudo.ordered.matrices,
       "pr_results"=procrust.results, "pr_matrices"=procrust.ordered.matrices, "swap"=SWAP, "lead"=lead, "scnd"=scnd)
}






#' Assessment of matrix similarity using several different types of metrics
#' Approaches include pearson and kendall correlation, classification error (Precision-recall curve),
#' categorization similarity (Cohen's Kappa metric), and Euclidian distance (frobenius norm)
#' All of these require columns to be mapped to each other, which can be done either by a greedy correlation search, or by a procrustes analysis
#' This function automatically does both
#' @param A
#' @param B
#' @param corr.type
#' @param full.procrust
#'
#' @return
#' @export
#'
#' @examples
compareModelMatricesComprehensive <- function(A,B, corr.type = "pearson")
{
  #First, align the matrices to have the same dimensions
  library(PRROC)
  prepped <- prepMatricesForAnalysis(A,B)
  lead <- prepped$lead;  scnd <- prepped$second; SWAP <- prepped$swap
  stopifnot(!SWAP)
  #Use my homeade correlation matching
  pseudo.procrust.results <- corrBasedScoring(lead, scnd, corr.method = "pearson")
  pseudo.ordered.matrices <- pseudo.procrust.results$corr_ordered_matrices

  #Use a published procrustes package too
  procrust.results <- procrustesBasedScoring(lead,scnd, SWAP)
  procrust.ordered.matrices <- procrust.results$procrust_ordered_matrices

  #Need to verify and check orders are consistent here, before we do downstream analysis
  #checkUVAlignment
  #All we need to update here- lead and scnd, x_pearson, x_kendall, x_rss,
  #  corr_based_stats <- with(comp.u[[n]],

  #repairUVAlignment <- function(v.dat, u.dat, original.fg.v, original.fg.u,original.uk.v,original.uk.u, source="procrust")

  #assess with greedy correlation based order
  by_pseudo_procrust_order <- classificationScoring(pseudo.ordered.matrices$lead,pseudo.ordered.matrices$scnd, SWAP)
  by_true_procrust_order <- classificationScoring(procrust.ordered.matrices$lead,
                                                  procrust.ordered.matrices$scnd, SWAP)
  names(by_true_procrust_order) <- gsub(names(by_true_procrust_order), pattern = "^", replacement = "procrust_")
  names(by_pseudo_procrust_order) <- gsub(names(by_pseudo_procrust_order), pattern = "^", replacement = "corr_based_")

  #Update results if all entries are 0:
  if((all(scnd == 0) & !all(lead == 0)) | (all(lead == 0) & !all(scnd == 0)))
  {
    message("An empty factorization appeared. Correlation scores being set to 0.")
    ret.list$procrustes.corr=0;ret.list$kendall.corr=0;ret.list$true.procrustes.corr=0
  }

  #Update the return data
  ret.list <- c(pseudo.procrust.results,procrust.results,by_pseudo_procrust_order,by_true_procrust_order )
  #How closely do our versions compare?
  ret.list$methods_mismatch <- compareProcrustesAndHomeade(pseudo.procrust.results$corr_based_data, procrust.results$procrustes_data,A,B,lead,scnd,DEBUG=FALSE)
  if(SWAP)
  {
    ret.list$order = "swapped";  ret.list$A = scnd;  ret.list[["B"]]=lead;
  }else
  {
    ret.list$order = "preserved";  ret.list$A = lead;  ret.list[["B"]]=scnd;
  }
  ret.list
}


#' Assessment of matrix similarity as clasfifiers
#' @param mat.dat
#' @param corr.type
#'
#' @return
#' @export
#'
#' @examples
compareMatricesAsClassifiersByMatrix <- function(mat.dat,A,B, corr.type = "pearson")
{
  #First, align the matrices to have the same dimensions
  library(PRROC)
  #lead <- mat.dat$lead;  scnd <- mat.dat$scnd; SWAP <- mat.dat$swap
  #  list("cb_results"=pseudo.procrust.results, "cb_matrices"=pseudo.ordered.matrices,
  #"pr_results"=procrust.results, "pr_matrices"=procrust.ordered.matrices, "swap"=SWAP, "lead"=lead, "scnd"=scnd)
  #assess with greedy correlation based order

  by_pseudo_procrust_order <- compareMatricesAsClassifiers(mat.dat$cb_matrices, "corr_based_")
  by_true_procrust_order <- compareMatricesAsClassifiers(mat.dat$pr_matrices, "procrust_")
  #by_pseudo_procrust_order <- classificationScoring(mat.dat$cb_matrices$lead, mat.dat$cb_matrices$scnd, SWAP)
  #by_true_procrust_order <- classificationScoring(mat.dat$pr_matrices$lead,
  #                                                mat.dat$pr_matrices$scnd, SWAP)
  #names(by_true_procrust_order) <- gsub(names(by_true_procrust_order), pattern = "^", replacement = "procrust_")
  #names(by_pseudo_procrust_order) <- gsub(names(by_pseudo_procrust_order), pattern = "^", replacement = "corr_based_")

  #Update the return data
  #too much going on htere, need to simplify
  ret.list <- c(mat.dat$pr_results,mat.dat$cb_results,
                by_pseudo_procrust_order,by_true_procrust_order )
  #ret.list <- list(by_pseudo_procrust_order,by_true_procrust_order )
  #Update results if all entries are 0:
  if((all(A == 0) & !all(B == 0)) | (all(A == 0) & !all(B == 0)))
  {
    message("An empty factorization appeared. Correlation scores being set to 0.")
    ret.list$procrustes.corr=0;ret.list$kendall.corr=0;ret.list$true.procrustes.corr=0
  }

  #How closely do our versions compare?
  ret.list$methods_mismatch <- compareProcrustesAndHomeade(mat.dat$cb_results$corr_based_data, mat.dat$pr_results$procrustes_data,A,B,DEBUG=FALSE)
  if(mat.dat$swap)
  {
    ret.list$order = "swapped";  ret.list$A = A;  ret.list[["B"]]=B;
  }else
  {
    ret.list$order = "preserved";  ret.list$A = A;  ret.list[["B"]]=B;
  }
  ret.list
}


#' Assessment of matrix similarity as clasfifiers
#' @param mat.dat
#' @param corr.type
#'
#' @return
#' @export
#'
#' @examples
compareMatricesAsClassifiersJoint <- function(mat.dat, corr.type = "pearson")
{
  #First, align the matrices to have the same dimensions
  library(PRROC)
  #first V
  by_corr_order.v <- compareMatricesAsClassifiers(list("lead"=mat.dat$cb_matrices$lead_v, "scnd"=mat.dat$cb_matrices$scnd_v, "swap"= mat.dat$cb_matrices$swap), "corr_based_")
  by_procrust_order.v <- compareMatricesAsClassifiers(list("lead"=mat.dat$pr_matrices$lead_v, "scnd"=mat.dat$pr_matrices$scnd_v, "swap"= mat.dat$cb_matrices$swap),
                                                    "procrust_")
  #then U
  by_corr_order.u <- compareMatricesAsClassifiers(list("lead"=mat.dat$cb_matrices$lead_u, "scnd"=mat.dat$cb_matrices$scnd_u, "swap"= mat.dat$cb_matrices$swap), "corr_based_")
  by_procrust_order.u <- compareMatricesAsClassifiers(list("lead"=mat.dat$pr_matrices$lead_u, "scnd"=mat.dat$pr_matrices$scnd_u, "swap"= mat.dat$cb_matrices$swap),
                                                      "procrust_")

  ret.list <- list("by_corr"=list("V"=by_corr_order.v, "U"=by_corr_order.u),
                   "by_procrust"=list("V"=by_procrust_order.v, "U"=by_procrust_order.u))

  #May wish to compare procrust and homeade here
  ret.list
}




###Reorganize above code to be more compartmentalized, transferable, whatever
#' Assessment of matrix similarity as clasfifiers
#' @param mat.dat - list containing lead, second and swap data to assess.
#' @param corr.type
#'
#' @return
#' @export
#'
#' @examples
compareMatricesAsClassifiers <- function(mat.dat, order_name_prefix)
{
  #First, align the matrices to have the same dimensions
  library(PRROC)
  lead <- mat.dat$lead;  scnd <- mat.dat$scnd; SWAP <- mat.dat$swap
  if(!is.null(SWAP)) {stopifnot(!SWAP)}else{SWAP <- FALSE}


  by_order <- classificationScoring(mat.dat$lead, mat.dat$scnd, SWAP)
  names(by_order) <- gsub(names(by_order), pattern = "^", replacement = order_name_prefix)

  #Update the return data
  ret.list <- by_order

  if(SWAP)
  {
    ret.list$order = "swapped";  ret.list$A = scnd;  ret.list[["B"]]=lead;
  }else
  {
    ret.list$order = "preserved";  ret.list$A = lead;  ret.list[["B"]]=scnd;
  }
  ret.list
}



#' Put the larger matrix first and fill empty columns with 0s
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
prepMatricesForAnalysis <- function(A,B)
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

#' Perform scoring on 2 matrices based on classification
#'
#' @param lead the alleged reference matrix (larger of the 2), ORDERED
#' @param scnd the alleged predicted matrix, ORDERED
#' @param SWAP if initial A and B correspond to lead and scnd, this is FALSE
#'
#' @return
#' @export
#'
#' @examples
classificationScoring <- function(lead, scnd, SWAP)
{
  #do a simple BINARY precision variance assessment
  prd <-getBinaryPRCurve(lead, scnd)

  #Do a multi-class assessment, looking at +/0/-
  kappa.score <- threeClassKappaComplete(lead, scnd)
  kappa.score.backwards <- threeClassKappaComplete(scnd, lead)
  if(!is.na(kappa.score$global$kappa)) #fails if its NA
  {
    stopifnot(kappa.score$global$kappa == kappa.score.backwards$global$kappa)
  }


  ret.list <- list("multi.class"=kappa.score)

  if(SWAP)
  {
    ret.list$order = "swapped";
    ret.list[["A.optimal"]] = scnd;  ret.list[["B.optimal"]]=lead
    ret.list[["caret.report.A.first"]]=prd$caret.report.bkwd; ret.list[["caret.report.B.first"]]=prd$caret.report.fwd
    ret.list[["metrics.A.first"]] = prd$metrics.bwkd;  ret.list[["metrics.B.first"]]=prd$metrics.fwd
    ret.list[["pr.A.first"]]=prd$pr.second; ret.list[["pr.B.first"]]=prd$pr.first
  }else
  {
    ret.list$order = "preserved";
    ret.list[["A.ordered"]] =lead;  ret.list[["B.ordered"]]=scnd
    ret.list[["caret.report.A.first"]]=prd$caret.report.fwd; ret.list[["caret.report.B.first"]]=prd$caret.report.bkwd
    ret.list[["metrics.A.first"]] =prd$metrics.fwd;  ret.list[["metrics.B.first"]]=prd$metrics.bwkd
    ret.list[["pr.A.first"]]=prd$pr.first; ret.list[["pr.B.first"]]=prd$pr.second
  }

  ret.list
}

#' Title
#'
#' @param lead alleged reference matrix, unordered
#' @param scnd alleged predicted matrix, unordered
#'
#' @return Results are with the transformation and turning, NOT WITH THE SPARSE VERSION OF THE MATRIX.
#' @export
#'
#' @examples
procrustesBasedScoring <- function(lead,scnd, DEBUG=FALSE)
{
  procestes.corr <- fullProcrustes(lead, scnd) #This is too slow.

    procrust.ordered.matrices <- list("lead"=lead,
                                    "scnd"=matrixSignsProduct(scnd[,procestes.corr$mapping_reorder],procestes.corr$mapping_signs))
  cor.by.procrust <- stackAndAssessCorr(procestes.corr$A.normalized, procestes.corr$B.transformed) #correlation with the raw adjusted matrices
  rss.by.procrust <- rrmse(procestes.corr$B.transformed,procestes.corr$A.normalized) #ORDER IS PRED, TRUE
  #rss.by.procrust <- rrmse(procrust.ordered.matrices$lead, procrust.ordered.matrices$scnd)
  corr.k.by.procrust <- stackAndAssessCorr(procestes.corr$A.normalized, procestes.corr$B.transformed, cor.type = "kendall")

  list("procrustes_data"=procestes.corr, "procrust_ordered_matrices" = procrust.ordered.matrices,
       "procrustes_pearson"=cor.by.procrust, "procrustes_rss"=rss.by.procrust, "procrustes_kendall"=corr.k.by.procrust)

}

#' Correlation and RSS scoring based on greedy correlation matching
#'
#' @param lead alleged reference matrix, unordered
#' @param scnd alleged predicted matrix, unordered
#' @param corr.method - which correlation to use to order factors, default is pearson
#'
#' @return
#' @export
#'
#' @examples
corrBasedScoring <- function(lead, scnd, corr.method = "pearson")
{
  #G et the match and the greedy correlation corresponding to it
  #Pseudo-procrustes is based on correlation, rather than euclidian distance
  aligned.corr <- pseudoProcrustes(lead, scnd,corr.type = corr.method) #TODO:check this is reasonably optimal.
  pseudo.procestes.corr <- aligned.corr$greedy.match$corr
  corr_ordered_matrices = list("lead"=aligned.corr$A, "scnd"=aligned.corr$B)
  #alternative: kendall:
  #aligned.corr.k <- pseudoProcrustes(lead, scnd,corr.type = "kendall") #TODO:check this is reasonably optimal. #this is much harder to calculate for U, weird still slow
  #^this version reorders by kendall
  #Below- retain the pearson-based order
  aligned.corr.k <- stackAndAssessCorr(aligned.corr$A, aligned.corr$B, cor.type = "kendall") #This step is slow for some reason....

  #RSS <- norm(aligned.corr$A - aligned.corr$B,  type = "F")
  RSS <- rrmse(aligned.corr$B,aligned.corr$A) #ORDER IS PRED, TRUE
  list("corr_based_data"=aligned.corr, "corr_ordered_matrices" = corr_ordered_matrices,
       "corr_based_pearson"= pseudo.procestes.corr, "corr_based_rss"=RSS, "corr_based_kendall"=aligned.corr.k)
}

#' Simple check to see if our correlation based match vs the prorustes based one yield the same alignment
#'
#' @param corr.alignment
#' @param procrustes.alignment
#' @param lead
#' @param scnd
#' @param DEBUG
#'
#' @return
#' @export
#'
#' @examples
compareProcrustesAndHomeade <- function(corr.alignment, procrustes.alignment,lead,scnd,DEBUG=TRUE)
{
  ret.dat <- FALSE
  A <- dropEmptyCols(lead)
  B <- dropEmptyCols(scnd)
  #Check- how closely does my alignment match with the procrustes version
  non.matching <- sum(corr.alignment$greedy.match$order.pred != procrustes.alignment$mapping_reorder)
  k.diff <- abs(ncol(A) - ncol(B)) #possible mismatches due to missing columns
  if(non.matching > k.diff)
  {
    message("In this test, our greedy approach yields a different ordering than true procrustes.")
    message("Review closely all possible outputs")
    ret.dat <- TRUE
  }
  if(DEBUG) #Look more closely at the correlation differences, don't need anymore
  {
    #Reorder according to procrustes
    pred.one <- as.matrix(scnd[,procestes.corr$mapping_reorder]) %*% diag(procestes.corr$mapping_signs);      lead.one <- lead
    cor.reorder <- stackAndAssessCorr(lead.one, pred.one)
    #Try the different possible settings
    if(cor.reorder > pseudo.procestes.corr)
    {
      message("Procrustes orientation yields a better correlation than ours.")
    }else if(cor.reorder < pseudo.procestes.corr){
      message("Our version outperforms procrustes analysis! No change to make")
    }else
    {
      message("The correlations are the same;, no change to make.")
    }
  }
  ret.dat
}

#Get a binary PRcurve by looking at the columsn of the lead and second maftrix
#Binarized such that 0/1
getBinaryPRCurve <- function(lead, scnd)
{
  #BAD ashton- don't fix code like this.
  #return(list("caret.report.bkwd"=NA,"caret.report.fwd"=NA,
  #            "metrics.bwkd" = NA, "metrics.fwd"=NA, "pr.second"=NA, "pr.first"=NA))
  library(caret)
  if(any(is.na(lead))|any(is.na(scnd)))
  {
    message("Found some NA, not good")
    print(lead)
    print(scnd)
  }
  lead.binarized <- lead
  scnd.binarized <- scnd
  lead.binarized[lead.binarized != 0] <- 1
  scnd.binarized[scnd.binarized != 0] <- 1
  predicted <- as.vector(scnd.binarized)
  actual <- as.vector(lead.binarized)
  #Using one first
  xtab <- table(predicted, actual)
  #print(xtab)
  #If entries are all ones, than skip this
  if(all(predicted == 1) | all(actual == 1) | nrow(xtab) < 2 | ncol(xtab) < 2 | all(predicted == 0) | all(actual == 0) )
  {
  #Special case to deal with
    if(all(predicted==1) & !all(actual ==1))
    {
      warning("Special case where all classified the same-- cannot calculate confusion matrix")
    }
    conf.forward <- NA;
    conf.backward <- NA
  }else
  {
    conf.forward <- caret::confusionMatrix(table(predicted, actual))
    conf.backward <- caret::confusionMatrix(table(actual, predicted))
  }

  #alternatively
  library(Metrics)
  forward.metrics <- list("recall"=Metrics::recall(actual = actual, predicted = predicted),
                          "precision" = Metrics::precision(actual = actual, predicted = predicted),
                          "f1"=   Metrics::f1(actual, predicted))

  reverse.metrics <- list("recall"=Metrics::recall(actual = predicted, predicted = actual),
                          "precision" = Metrics::precision(actual = predicted, predicted = actual),
                          "f1"=   Metrics::f1(predicted, actual))

  #PR curve:
  if(!all(scnd == 0))
  {
    scnd.scores <- scnd^2/max(scnd^2)
  }
  else
  {
    empt.list <- list("type"="PR", "auc.integral"=NA,"auc.davis.goadrich"=NA)
    return(list("caret.report.bkwd"=conf.backward,"caret.report.fwd"=conf.forward,
                "metrics.bwkd" = reverse.metrics, "metrics.fwd"=forward.metrics, "pr.second"=empt.list, "pr.first"=empt.list))
  }


  #scores of all the trues.
  positives <- scnd.scores[lead.binarized != 0]
  negatives <- scnd.scores[lead.binarized == 0]
  curve.lead.first <- pr.curve(scores.class0 = positives, scores.class1 = negatives, curve = TRUE)


  lead.scores <- lead^2/max(lead^2)
  if(all(lead ==0 ))
  {
    lead.scores <- lead
  }
  #scores of all the trues.
  if(all(lead.scores == 0))
  {
    message("All scores are 0. Randomly making one drawm from a uniform to allow for non-error")
    warning("All scores are 0. Randomly making one drawm from a uniform to allow for non-error")
    lead.scores[sample(1:length(lead.scores), size = 1)] <- runif(1)

  }
  positives <- lead.scores[scnd.binarized != 0]
  negatives <- lead.scores[scnd.binarized == 0]
  curve.scnd.first <- pr.curve(scores.class0 = positives, scores.class1 = negatives, curve = TRUE)


  return(list("caret.report.bkwd"=conf.backward,"caret.report.fwd"=conf.forward,
              "metrics.bwkd" = reverse.metrics, "metrics.fwd"=forward.metrics, "pr.second"=curve.scnd.first, "pr.first"=curve.lead.first))

}

splitThreeClasses <- function(v)
{
  lead.classes <- v
  lead.classes[lead.classes > 0] <- 1
  lead.classes[lead.classes < 0] <- 2
  lead.classes
}
#TODO: verify this is same regardless of order (it should be, since looking at similarity)
threeClassKappa <- function(lead, scnd, by_factor = FALSE)
{
  lead.classes <- splitThreeClasses(lead)
  scnd.classes <- splitThreeClasses(scnd)
  lead.vector <- as.vector(lead.classes) #(stacks columns)
  scnd.vector <- as.vector(scnd.classes)
  global.scores <- suppressMessages(psych::cohen.kappa(data.frame(factor(lead.vector), factor(scnd.vector))))
  #message("Now by factor...")
  #Do it by column too
  #Be warned- this results in very sparse setups and so can yield error pars that are invalid.
  by.col = NA
  if(by_factor)
  {
    #BEWARE- suppressing warnings in this part.
    by.col <- lapply(1:ncol(lead.classes), function(i) suppressWarnings(psych::cohen.kappa(data.frame(lead.classes[,i], scnd.classes[,i]))))
  }

  return(list("global"=global.scores, "by.factor"=by.col))
}
##Its possible that the signs are inverted, which would result in poor scoring. Pick the orientation that leads to best
threeClassKappaComplete <- function(lead, scnd,...)
{
  pos <- threeClassKappa(lead, scnd, ...)
  neg <- threeClassKappa(lead, -1 * scnd, ...)
  if(neg$global$kappa > pos$global$kappa)
  {
    #message("Inverting sign on kappa score.")
    return(neg)
  }else
  {
    return(pos)
  }
}


#Do we run across all the factorizations, or do we just specify the falls?
#Do across all, and get the stats in an output file.
#Or maybe, even better, write some functions to call from a markdown notebook where you actualy make the figures, eh?

#read in all the finngen and ukbb results into lists

loadGLEANERResults <- function(dir, data_source)
{
  finngen.list <- list()
  for(f in list.files(path = dir, pattern = paste0("*", data_source)))
  {
    load(paste0(dir, f))
    name = gsub(x = f, replacement = "", pattern = data_source)
    if(data_source == "_final_dat.RData")
    {
      finngen.list[[name]] <- ret
    }else if(data_source == "_bic_dat.RData")
    {
      finngen.list[[name]] <- bic.dat
      finngen.list[[name]]$V <- finngen.list[[name]]$optimal.v
      finngen.list[[name]]$U <- finngen.list[[name]]$optimal.u
    }else
    {
      message("Invalid dataset. Returning")
      return(NA)
    }

  }
  return(finngen.list)
}

loadGLEANERResultsLists <- function(dir_finngen, dir_ukbb, dat_source="_final_dat.RData")
{
  finngen.list <- loadGLEANERResults(dir_finngen,dat_source)
  ukbb.list <- loadGLEANERResults(dir_ukbb,dat_source)
  return(list("finngen" = finngen.list, "ukbb"=ukbb.list))
}

#Perform the model comparison on a given matrix pair from finngen.list and ukbb.list
#n : the arg name to look at
#BAD CODING PRACTICE- we add directly to the table in this function and then update it accordingly.
#compareFactPrecision(n, finngen.list, ukbb.list, comp.v,comp.u,global.df.V,global.df.U,full.df.V,full.df.U)
#Change on Aug 16- remove all SWAP. Never swap, always reference by finngen (first one)
#Classifiers try both ways, but otherwise its always with respect to the first one
compareFactPrecision <- function(n, finngen.list, ukbb.list, comp.v,comp.u, global.df.V,global.df.U, overall.df.V, overall.df.U,global.joint, check.realign=TRUE)
{
  #Prepare data for analysis
  fg.v <- as.matrix(finngen.list[[n]]$V)
  uk.v <- as.matrix(ukbb.list[[n]]$V)
  fg.u <- as.matrix(finngen.list[[n]]$U)
  uk.u <- as.matrix(ukbb.list[[n]]$U)
  #Make sure the matrices have aligned dimensions before performing all the tests:
  v.dat <- prepMatricesForAnalysis(fg.v,uk.v);
  u.dat <- prepMatricesForAnalysis(fg.u,uk.u)
  stopifnot(!v.dat$swap)
  fg.v <- v.dat$lead; uk.v <- v.dat$second
  fg.u <- u.dat$lead; uk.u <- u.dat$second

  #Evaluation and alignment of matrices jointly
  forced.pairing <- forcedPairedEvaluation(fg.u,fg.v,uk.u,uk.v)
  joined.classification.data <- compareMatricesAsClassifiersJoint(forced.pairing)
  global.joint <- tabularizePairedData(forced.pairing$numeric_dat, joined.classification.data,global.joint)
  #Evaluation and alignment of matrices separately
  aligned.v <- alignMatrices(fg.v, uk.v)
  aligned.u <- alignMatrices(fg.u, uk.u)
  if(check.realign)
  {
    align.list <- ensureUVAlign(aligned.u,aligned.v,
                                aligned.u$lead,aligned.u$scnd,aligned.v$lead,aligned.v$scnd)
    aligned.v <- align.list$V; aligned.u <- align.list$U
  }else
  {
    align.list <- list(); align.list$alignment <- NA
  }

  #Verify the update that was needed was made...

  comp.v[[n]] <- compareMatricesAsClassifiersByMatrix(aligned.v, fg.v, uk.v)
  comp.u[[n]] <- compareMatricesAsClassifiersByMatrix(aligned.u, fg.u, uk.u)
  #TODO- combine all the data for easy parsing downhere....
  #comp.v[[n]] <- compareModelMatricesComprehensive(fg.v,uk.v)
  #  comp.u[[n]] <- compareModelMatricesComprehensive(fg.u,uk.u)
  #statistics overall
  cols.corr_based = c("cb_recall", "cb_precision", "cb_pearson", "cb_kendall", "cb_rss", "cb_auc", "cb_kappa")
  corr_based_stats <- with(comp.v[[n]],
                           c(corr_based_metrics.A.first$recall,corr_based_metrics.A.first$precision,
                             corr_based_pearson,corr_based_kendall,corr_based_rss,
                             corr_based_pr.A.first$auc.integral,corr_based_multi.class$global$kappa))

  cols.procrustes = c("pr_recall", "pr_precision", "pr_pearson", "pr_kendall", "pr_rss", "pr_auc", "pr_kappa")
  procrust_based_stats <- with(comp.v[[n]],
                           c(procrust_metrics.A.first$recall,procrust_metrics.A.first$precision,
                             procrustes_pearson,procrustes_kendall,procrustes_rss,
                             procrust_pr.A.first$auc.integral, procrust_multi.class$global$kappa))
  global.df.V <- cbind(global.df.V, c(paste0(ncol(fg.v), "-",ncol(uk.v)), corr_based_stats, procrust_based_stats))
  rownames(global.df.V) <- c("Fg_K-UK_K",cols.corr_based,cols.procrustes)

  #Statistics by factor (DEPRECATED)
  #overall.df.V <- rbind(overall.df.V,data.frame("study" = n, "Corr"=comp.v[[n]]$factor.procrust.corr, "Kappas"=sapply(comp.v[[n]]$multi.class$by.factor, function(x) x$kappa)))
  overall.df.V <- NA

  corr_based_stats <- with(comp.u[[n]],
                           c(corr_based_metrics.A.first$recall,corr_based_metrics.A.first$precision,
                             corr_based_pearson,corr_based_kendall,corr_based_rss,
                             corr_based_pr.A.first$auc.integral,corr_based_multi.class$global$kappa))

  procrust_based_stats <- with(comp.u[[n]],
                               c(procrust_metrics.A.first$recall,procrust_metrics.A.first$precision,
                                 procrustes_pearson,procrustes_kendall,procrustes_rss,
                                 procrust_pr.A.first$auc.integral, procrust_multi.class$global$kappa))

  global.df.U <- cbind(global.df.U, c(paste0(ncol(as.matrix(finngen.list[[n]]$V)), "-",ncol(as.matrix(ukbb.list[[n]]$V))), corr_based_stats, procrust_based_stats))
  rownames(global.df.U) <- c("Fg_K-UK_K",cols.corr_based,cols.procrustes)
    #Deprecated
  #overall.df.U <- rbind(overall.df.U,data.frame("study" = n,
  #                                              "Corr"=comp.u[[n]]$factor.procrust.corr,
  #                                              "Kappas"=sapply(comp.u[[n]]$multi.class$by.factor, function(x) x$kappa)))
  overall.df.U <- NA

  return(list(comp.v, global.df.V, comp.u, global.df.U,overall.df.V,overall.df.U, align.list$alignment,global.joint)) #last entry-alignment status.
}



#Paired stats, structure is a bit different
tabularizePairedData <- function(numeric.dat, classification.dat,global.joint)
{
  #rrmse is caled in the case of cb.
  class.cols <- c("joint_##_recall_@@", "joint_##_precision_@@", "joint_##_auc_@@", "joint_##_kappa_@@", "joint_##_pearson_@@","joint_##_kendall_@@","joint_##_rrmse_@@")
  cols.corr_based.v = gsub(class.cols, pattern = "##", replacement = "cb") %>% gsub(., pattern="@@", replacement = "v")
  corr_based_stats_joint.v <- c(with(classification.dat$by_corr$V,
                                     c(corr_based_metrics.A.first$recall,corr_based_metrics.A.first$precision,
                                       corr_based_pr.A.first$auc.integral,corr_based_multi.class$global$kappa)),
                                numeric.dat$cb_corr_v, numeric.dat$cb_kendall_v,numeric.dat$cb_rrmse_v )

  cols.pr_based.v = gsub(class.cols, pattern = "##", replacement = "pr") %>% gsub(., pattern="@@", replacement = "v")
  procrust_based_stats_joint.v <- c(with(classification.dat$by_procrust$V,
                                         c(procrust_metrics.A.first$recall,procrust_metrics.A.first$precision,
                                           procrust_pr.A.first$auc.integral,procrust_multi.class$global$kappa)),
                                    numeric.dat$pearson_v_paired, numeric.dat$kendall_v_paired,numeric.dat$rmse_v_paired)

  cols.corr_based.u = gsub(class.cols, pattern = "##", replacement = "cb") %>% gsub(., pattern="@@", replacement = "u")
  corr_based_stats_joint.u <- c(with(classification.dat$by_corr$U,
                                     c(corr_based_metrics.A.first$recall,corr_based_metrics.A.first$precision,
                                       corr_based_pr.A.first$auc.integral,corr_based_multi.class$global$kappa)),
                                numeric.dat$cb_corr_u, numeric.dat$cb_kendall_u,numeric.dat$cb_rrmse_u )

  cols.pr_based.u = gsub(class.cols, pattern = "##", replacement = "pr") %>% gsub(., pattern="@@", replacement = "u")
  procrust_based_stats_joint.u <- c(with(classification.dat$by_procrust$V,
                                         c(procrust_metrics.A.first$recall,procrust_metrics.A.first$precision,
                                           procrust_pr.A.first$auc.integral,procrust_multi.class$global$kappa)),
                                    numeric.dat$pearson_u_paired, numeric.dat$kendall_u_paired,numeric.dat$rmse_u_paired)

  all.names <- c(cols.corr_based.v,cols.pr_based.v,cols.corr_based.u,cols.pr_based.u)
  dat <- c(corr_based_stats_joint.v,procrust_based_stats_joint.v,corr_based_stats_joint.u,procrust_based_stats_joint.u)
  global.joint <- cbind(global.joint, dat)
  rownames(global.joint) <- all.names
  return(global.joint)
}









checkUVAlignment <- function(v.dat, u.dat, fg.v, fg.u,uk.v,uk.u, source="procrust")
{
  if(source == "procrust")
  {
    mismatches <- sum(v.dat$procrustes_data$mapping_reorder != u.dat$procrustes_data$mapping_reorder)
  }
  if(source == "corr_based")
  {
    mismatches <- sum(v.dat$corr_based_data$greedy.match$order.pred != u.dat$corr_based_data$greedy.match$order.pred)
  }

  k.off <- abs(ncol(fg.v) - ncol(uk.v)) #possible mismatches with 0 columns
  if(mismatches > k.off)
  {
    warning("U and V weren't aligned in the same fashion, you have been warned.")
    if(source == "corr_based")
    {
      test.fix <- repairUVAlignment(v.dat, u.dat, fg.v, fg.u,uk.v,uk.u, source=source)
    }

    return(TRUE)
  }
  return(FALSE)

}


repairUVAlignment <- function(v.dat, u.dat, original.fg.v, original.fg.u,original.uk.v,original.uk.u, source="procrust")
{
  if(source == "procrust")
  {
    mismatches <- sum(v.dat$procrustes_data$mapping_reorder != u.dat$procrustes_data$mapping_reorder)
  }
  if(source == "corr_based")
  {
    mismatches <- sum(v.dat$corr_based_data$greedy.match$order.pred != u.dat$corr_based_data$greedy.match$order.pred)
  }
  k.off <- abs(ncol(original.fg.v) - ncol(original.uk.v))
  if(mismatches > k.off)
  {
    #an easy option- if the number of mismatched < 5, just use the other version....
    #update to have 0s in it, shoot.
    ret <- fillWithZeros(original.fg.v, original.uk.v)
    original.uk.v <- ret$fp
    ret <- fillWithZeros(original.fg.u, original.uk.u)
    original.uk.u <- ret$fp

    if(mismatches <= 5)
    {
      #first, identify the mismatched columns (based on index in updated matrix)
      mismatch.v <- which((v.dat$corr_based_data$greedy.match$order.pred != u.dat$corr_based_data$greedy.match$order.pred)) #index of mismatches in new space and in reference (fg)
      v.cols.off <- v.dat$corr_based_data$greedy.match$order.pred[mismatch.v] #get the columns #s in ORIGINAL space
      mismatch.u <- which((u.dat$corr_based_data$greedy.match$order.pred !=v.dat$corr_based_data$greedy.match$order.pred)) #index of mismatches in new space and in reference (U)
      u.cols.off <- u.dat$corr_based_data$greedy.match$order.pred[mismatch.u] #original column #s in U
      stopifnot(all(u.cols.off %in% v.cols.off)); stopifnot(length(u.cols.off) == length(v.cols.off)) #The column numbers should be the same in both cases, just ordered differently
      #extract those same columns from U and V and make sure they are in the right order
      sorted.off <- sort(u.cols.off) #Sort them so they are order consistent; indices from the ORIGINAL,unsorted matrix that correspond to the sorted matrix

      #Extract just the columns to re-score
      fg.mismatch.v <- original.fg.v[,mismatch.u]
      uk.mismatch.v <- original.uk.v[,sorted.off]
      fg.mismatch.u <- original.fg.u[,mismatch.u]
      uk.mismatch.u <- original.uk.u[,sorted.off]

      #Extract the columns that are already in a satisfactory order
      fg.match.v <- v.dat$corr_based_data$A[,-mismatch.v]
      uk.match.v <- v.dat$corr_based_data$B[,-mismatch.v]
      fg.match.u <- u.dat$corr_based_data$A[,-mismatch.u]
      uk.match.u <- u.dat$corr_based_data$B[,-mismatch.u]

      #Find the order that maximizes R2
      fixed.order <- evaluate_error(fg.mismatch.u, fg.mismatch.v, uk.mismatch.u, uk.mismatch.v)

      #Build the new matrices, with order
      uk.match.u <- cbind(uk.match.u,fixed.order$reordered_U)
      uk.match.v <- cbind(uk.match.v,fixed.order$reordered_V)

      fg.match.u <- cbind(fg.match.u,fg.mismatch.u)
      fg.match.v <- cbind(fg.match.v,fg.mismatch.v )

      #Get the new correlation:
      #update the order list to reflect the unification

      #Check that the order is now consistent.

      return(list("fg_v"=fg.match.v, "fg_u"=fg.match.u, "uk_v"=uk.match.v, "uk_u"=uk.match.u, "v_corr"=stackAndAssessCorr(fg.match.v,uk.match.v),"u_corr"= stackAndAssessCorr(fg.match.u,uk.match.u)))

    }else
    {
      message("not implemented yet....")
      reorder.u <- getNewOrder(u.dat)
      reorder.v <- getNewOrder(v.dat)
      u.order.rss <- distanceByOrder(original.fg.u, original.fg.v, original.uk.u, original.uk.v, reorder.u)
      v.order.rss <- distanceByOrder(original.fg.u, original.fg.v, original.uk.u, original.uk.v, reorder.v)
      if(u.order < v.order)
      {
        v.dat <- updateProcrustOrder(v.dat,reorder.u )
      }else
      {
        u.dat <- updateProcrustOrder(u.dat,reorder.v)
      }
      return(list(v.dat, u.dat))
    }
  }
}

#A:lead
#B:second
ensureUVAlign <- function(u.dat, v.dat,A_u,B_u, A_v, B_v)
{
  #If alignemtn needs doing....
  alignment.repairs <- detectAlignmentIssues(u.dat, v.dat,A_v, B_v) #For this, we need A,B with 0 columns dropped
  result = switch(
    alignment.repairs,
    "none"= list("U" = u.dat, "V" = v.dat),
    "procrustes"= repairUVAlignmentProcrustes(u.dat, v.dat,A_u,B_u, A_v, B_v),
    "corr_based"= repairUVAlignmentCorrBased(u.dat, v.dat,A_u,B_u, A_v, B_v),
    "both"= repairBoth(u.dat, v.dat,A_u,B_u, A_v, B_v)
  )
c(result, list("alignment"=alignment.repairs))
}

#' Repair alignment both in procrusts and correlation based famework. Baby helper function
#'
#' @param u.dat
#' @param v.dat
#' @param A_u
#' @param B_u
#' @param A_v
#' @param B_v
#'
#' @return
#' @export
repairBoth <- function(u.dat, v.dat,A_u,B_u, A_v, B_v)
{
  pf <- repairUVAlignmentProcrustes(u.dat, v.dat,A_u,B_u, A_v, B_v)
  repairUVAlignmentCorrBased(pf$U, pf$V,A_u,B_u, A_v, B_v)
}


#' Determine if we need to align matrices at all.
#'
#' @param v.dat
#' @param u.dat
#' @param original.fg.v
#' @param original.uk.v
#'
#' @return
#' @export
#'
#' @examples
detectAlignmentIssues <- function(v.dat, u.dat, original.fg.v, original.uk.v)
{
  original.fg.v <- dropEmptyCols(original.fg.v)
  original.uk.v <- dropEmptyCols(original.uk.v)
  mismatches.procrust <- sum(v.dat$pr_results$procrustes_data$mapping_reorder != u.dat$pr_results$procrustes_data$mapping_reorder)
  mismatches.corrbased <- sum(v.dat$cb_results$corr_based_data$greedy.match$order.pred != u.dat$cb_results$corr_based_data$greedy.match$order.pred)
  k.off <- abs(ncol(original.fg.v) - ncol(original.uk.v))
  if(mismatches.procrust > k.off & mismatches.corrbased > k.off)
  {
    return("both")
  }else if(mismatches.procrust > k.off)
  {
    return("procrustes")
  }else if(mismatches.corrbased > k.off)
  {
    return("corr_based")
  }else
  {
    return("none")
  }

}

dropEmptyCols <- function(X)
{
  non_empty = which(apply(X, 2, function(x) sum(x!=0)) > 0)
  as.matrix(X[,non_empty])
}
repairReorder <- function(u.dat,v.dat,A_u, A_v, B_u, B_v,version)
{
  reorder.u <- getNewOrder(u.dat, version)
  reorder.v <- getNewOrder(v.dat, version)
  #trueU, trueV, predU,predV,new.order
  u.order.rss <- distanceByOrder(A_u, A_v, B_u, B_v, reorder.u)
  v.order.rss <- distanceByOrder(A_u, A_v, B_u, B_v, reorder.v)

  u.order.r2 <- r2ByOrder(A_u, A_v, B_u, B_v, reorder.u)
  v.order.r2 <- r2ByOrder(A_u, A_v, B_u, B_v, reorder.v)
  return(list("u.rss" = u.order.rss, "v.rss" = v.order.rss, "u.r2" = u.order.r2,"v.r2"=v.order.r2, "reorder.u"=reorder.u, "reorder.v" = reorder.v))
}
repairUVAlignmentCorrBased<- function(u.dat, v.dat,A_u,B_u, A_v, B_v)
{
  #message("NOTE: need to realign U and V based on correlation here")
  r <- repairReorder(u.dat, v.dat,A_u,A_v, B_u, B_v, "corr_based")
  if(r$u.r2 > r$v.r2)
  {
    v.dat <- updateCBOrder(v.dat,B_v,r$reorder.u )
  }else
  {
    u.dat <- updateCBOrder(u.dat,B_u,r$reorder.v)
  }

  return(list("V" = v.dat, "U"=u.dat))
}


repairUVAlignmentProcrustes <- function(u.dat, v.dat,A_u,B_u, A_v, B_v)
{
  message("NOTE: need to realign U and V based on procrustes here")
  r <- repairReorder(u.dat, v.dat,A_u,A_v, B_u, B_v, "procrustes")
  #Note- "with" needs to return something.
 if(r$u.rss < r$v.rss)
       {
        v.dat <- updateProcrustOrder(v.dat, B_v, r$reorder.u )
      }else {
        u.dat <- updateProcrustOrder(u.dat, B_u, r$reorder.v)
      }
  return(list("V" = v.dat, "U"=u.dat))
}

 getNewOrder <- function(mat.dat,source)
 {
   #return the original matrix order from the object
   if(source == "procrustes")
   {
     list("order" = mat.dat$pr_results$procrustes_data$mapping_reorder,
          "signs" = mat.dat$pr_results$procrustes_data$mapping_signs)
     #v.dat$pr_results$procrustes_data$mapping_reorder
   }else
   {
     list("order" =mat.dat$cb_results$corr_based_data$greedy.match$order.pred,
          "signs" =mat.dat$cb_results$corr_based_data$greedy.match$signs)
   }

 }

 distanceByOrder <- function(trueU, trueV, predU,predV,new.order)
 {
   rrmse(matrixSignsProduct(predU[,new.order$order], new.order$signs),trueU) +
     rrmse(matrixSignsProduct(predV[,new.order$order], new.order$signs),trueV)
 }

 r2ByOrder <- function(trueU, trueV, predU,predV,new.order)
 {
   stackAndAssessCorr(trueU,matrixSignsProduct(predU[,new.order$order], new.order$signs))^2 +
     stackAndAssessCorr(trueV,matrixSignsProduct(predV[,new.order$order], new.order$signs))^2
 }

 #All we need to update here- lead and scnd, x_pearson, x_kendall, x_rss,
 #  library(PRROC)
 #prepped <- prepMatricesForAnalysis(A,B)
 #lead <- prepped$lead;  scnd <- prepped$second; SWAP <- prepped$swap

 #Use my homeade correlation matching
 #pseudo.procrust.results <- corrBasedScoring(lead, scnd, corr.method = "pearson")
 #pseudo.ordered.matrices <- pseudo.procrust.results$corr_ordered_matrices

 #Use a published procrustes package too
 #procrust.results <- procrustesBasedScoring(lead,scnd)
 #procrust.ordered.matrices <- procrust.results$procrust_ordered_matrices
 ##HERE HERE HERE
 updateProcrustOrder <- function(obj, base.matrix, new.order)
 {
   #  #  list("cb_results"=pseudo.procrust.results, "cb_matrices"=pseudo.ordered.matrices,
   #"pr_results"=procrust.results, "pr_matrices"=procrust.ordered.matrices, "swap"=SWAP)
   update <- obj
   #update matrices
   update$pr_matrices$scnd <- matrixSignsProduct(base.matrix[,new.order$order], new.order$signs)

   #update order
   update$pr_results$procrustes_data$mapping_reorder <- new.order$order
   update$pr_results$procrustes_data$mapping_signs <- new.order$signs
   #update corr (rho and pearson)
   update$pr_results$procrustes_pearson <- stackAndAssessCorr(update$pr_matrices$lead,update$pr_matrices$scnd,cor.type = "pearson")
   update$pr_results$procrustes_kendall <- stackAndAssessCorr(update$pr_matrices$lead,update$pr_matrices$scnd,cor.type = "kendall")
   #update RSS
   update$pr_results$procrustes_rss <- rrmse(update$pr_matrices$scnd,update$pr_matrices$lead) #Order is PRED, TRUE
   return(update)
 }


 updateCBOrder <- function(obj, base.matrix, new.order)
 {
   update <- obj
   #update matrices
   update$cb_matrices$scnd <- matrixSignsProduct(base.matrix[,new.order$order], new.order$signs)
   #corr.alignment$greedy.match$order.pred != procrustes.alignment$mapping_reorder
   #update order
   update$cb_results$corr_based_data$greedy.match$order.pred <- new.order$order
   update$cb_results$corr_based_data$greedy.match$signs <- new.order$signs
   #update corr (rho and pearson)
   update$cb_results$corr_based_pearson <- stackAndAssessCorr(update$cb_matrices$lead,update$cb_matrices$scnd,cor.type = "pearson")
   update$cb_results$corr_based_kendall <- stackAndAssessCorr(update$cb_matrices$lead,update$cb_matrices$scnd,cor.type = "kendall")
   #update RSS
   update$cb_results$corr_based_rss <- rrmse(update$cb_matrices$scnd,update$cb_matrices$lead) #ORDER IS PRED, tRUE
   return(update)
 }
##
#Pass in the lists of runs, and perform the comparison tests for each, and write the results out in an easy-to-read table.
#Note that this is not a fast function, the comparison step is kinda slow.
#Returns table with data for Recall, precision, correlation, K counts, and AUPRC
#makeTableOfComparisons(no.shrink.results$finngen, no.shrink.results$ukbb)
#Re-factor this code: give global scores as well as distributional scores across factors
#direct.gamma.comparisons <- makeTableOfComparisons(std_block.results$finngen, std_block.results$ukbb)
makeTableOfComparisons <- function(finngen.list, ukbb.list, filter_term = "BIC", by_order = FALSE,...)
{
  both <- intersect(names(ukbb.list), names(finngen.list))
  if(by_order)
  {
    both = 1:length(ukbb.list)
  }
  #drop BIC-sklearn_K-MAX_COVAR
  #both <- both[-which(both == "BIC-sklearn_K-MAX_COVAR")]
  comp.v <- list()
  comp.u <- list()
  comp.xhat <- list()
  comp.mimsatched <- c()
  #These get the global statistics, summarizing the entire factorization
  global.df.V <- NULL
  global.df.U <- NULL
  global.df.joint <- NULL
  #These get the statistics by factor, allowing us to make boxplots
  full.df.V <- NULL
  full.df.U <- NULL
  full.df.beta <- NULL
  forced.pairing.dat <- NULL
  missing.dat <- c()
  both.caught <- c()
  for(n in both)
  {
    if(!(grepl(pattern=filter_term, x= n)))
    {
      next
    }
    print(n)
    both.caught <- c(both.caught, n)
    if(length(ukbb.list[[n]]) == 0 | (length(finngen.list[[n]]) == 0))
    {
      missing.dat <- c(missing.dat, n)
      global.df.V <- cbind(global.df.V, rep(NA,nrow(global.df.V)))
      global.df.U <- cbind(global.df.U, rep(NA,nrow(global.df.U)))
      global.df.joint <- cbind(global.df.joint, rep(NA,nrow(global.df.joint)))
      comp.xhat<- c(comp.xhat,NA)
      next
    }
    nfd <- compareFactPrecision(n, finngen.list, ukbb.list, comp.v,comp.u,global.df.V,global.df.U,full.df.V,full.df.U,global.df.joint,...)
    comp.v<-nfd[[1]]; global.df.V<-nfd[[2]]; comp.u<-nfd[[3]]; global.df.U<-nfd[[4]];full.df.V <- nfd[[5]]; full.df.U<- nfd[[6]];
    #forced.pairing.dat <- rbind(forced.pairing.dat, unlist(nfd[[8]]$numeric_dat))
    global.df.joint <- nfd[[8]]
    #Comparison of X directly
    #Make sure the orders are the same.
    stopifnot(all(finngen.list[[n]]$trait.names == ukbb.list[[n]]$trait.names))
    stopifnot(all(finngen.list[[n]]$snp.ids == ukbb.list[[n]]$snp.ids))

    beta.dat <- compareEstimatedEffectSizes(finngen.list[[n]]$U, ukbb.list[[n]]$U,finngen.list[[n]]$V, ukbb.list[[n]]$V)
    comp.xhat<- c(comp.xhat, beta.dat$rrmse_beta); full.df.beta <- rbind(full.df.beta, beta.dat$classification_beta)

    comp.mimsatched <-  c(comp.mimsatched, nfd[[7]])
  }

  missing_i <- which(both %in% missing.dat)
  #rownames(global.df.V) <- c("recall", "precision", "correlation", "Fg_K-UK_K", "PRC", "GlobalKappa", "Kendall_correlation", "Procrustes_pearson", "Euclidian_dist")
  colnames(global.df.V)<- both.caught#[-missing_i]
  colnames(global.df.joint) <- both.caught
  #rownames(global.df.U) <- c("recall", "precision", "correlation", "Fg_K-UK_K", "PRC","GlobalKappa","Kendall_correlation", "Procrustes_pearson", "Euclidian_dist")
  colnames(global.df.U)<- both.caught#[-missing_i]

  colnames(full.df.beta) <- c("recall", "precision", "auc", "kappa")

  pr.dat.v <- data.frame(t(global.df.V)) %>% tibble::rownames_to_column("rn") %>%
    separate(rn, into = c("BIC", "Factors"), sep = "_K") %>%
    mutate("BIC" = gsub(replacement = "", pattern = "BIC-", x = BIC),"Factors" = gsub(replacement = "", pattern = "-", x = Factors)) %>%
    mutate_at(colnames(.)[-c(1,2,3)], as.numeric)
  #pr.dat.v$recall <- as.numeric(pr.dat.v$recall)
  #pr.dat.v$precision <- as.numeric(pr.dat.v$precision)
  #pr.dat.v$correlation <- as.numeric(pr.dat.v$correlation)
  #pr.dat.v$PRC <- as.numeric(pr.dat.v$PRC)

   pr.dat.u <- data.frame(t(global.df.U)) %>% tibble::rownames_to_column("rn") %>%
    separate(rn, into = c("BIC", "Factors"), sep = "_K") %>%
    mutate("BIC" = gsub(replacement = "", pattern = "BIC-", x = BIC),"Factors" = gsub(replacement = "", pattern = "-", x = Factors))  %>%
     mutate_at(colnames(.)[-c(1,2,3)], as.numeric)

   joint.df.return <- data.frame(t(global.df.joint)) %>% tibble::rownames_to_column("rn") %>%
     separate(rn, into = c("BIC", "Factors"), sep = "_K") %>%
     mutate("BIC" = gsub(replacement = "", pattern = "BIC-", x = BIC),"Factors" = gsub(replacement = "", pattern = "-", x = Factors))  %>%
     mutate_at(colnames(.)[-c(1,2,3)], as.numeric)
   #forced.pairing.dat <- cbind(pr.dat.u %>% select(BIC, Factors),forced.pairing.dat )

   #%>% print()
  #pr.dat.u$recall <- as.numeric(pr.dat.u$recall)
  #pr.dat.u$precision <- as.numeric(pr.dat.u$precision)
  #pr.dat.u$correlation <- as.numeric(pr.dat.u$correlation)
  #pr.dat.u$PRC <- as.numeric(pr.dat.u$PRC)

  #Directly compare the output X matrices. This may be a better metric
    #Add some more comprehensive data tables too
  return(list("V"=pr.dat.v, "U"=pr.dat.u, "X_hat.norms"=comp.xhat,"U-V-misaligned"=comp.mimsatched,
              "full.V" =full.df.V, "full.U" = full.df.U, "matrix_paired_results"=joint.df.return, "df_beta"=full.df.beta))

}


compareEstimatedEffectSizes <- function(fg_u, uk_u, fg_v, uk_v)
{
  Xhat.fg <- fg_u %*% t(fg_v)
  Xhat.uk <-uk_u %*% t(uk_v)

  if(any(dim(Xhat.uk) != dim(Xhat.fg)))
  {
    if(any(is.na(Xhat.fg)))
    {
      message("fg matrix is empty")
      Xhat.fg <- matrix(0,nrow=nrow(Xhat.uk), ncol=ncol(Xhat.uk))
    }
    if(any(is.na(Xhat.uk)))
    {
      message("uk matrix is empty")
      Xhat.uk <- matrix(0,nrow=nrow(Xhat.fg), ncol=ncol(Xhat.fg))
    }
  }
  #updating to be the root mean squared error:
  #comp.xhat<- c(comp.xhat, sqrt((norm(Xhat.fg - Xhat.uk, "F")^2)/ (ncol(Xhat.fg) * nrow(Xhat.fg))))
  comp.xhat<-  rrmse(Xhat.uk,Xhat.fg) #ORDER IS PRED, TRUE

  #First, align the matrices to have the same dimensions
  library(PRROC)
  simple.class.scoring <- classificationScoring(Xhat.fg, Xhat.uk, FALSE)

  #cols.beta_based = c("cb_recall", "cb_precision", "cb_pearson", "cb_kendall", "cb_rss", "cb_auc", "cb_kappa")
  beta_based_stats <- with(simple.class.scoring,
                           c(metrics.A.first$recall,metrics.A.first$precision,
                             pr.A.first$auc.integral,multi.class$global$kappa))

 return(list("rrmse_beta" = comp.xhat, "classification_beta"=beta_based_stats))

}




#' Compare all the covariance-adjusted outputs pairwise
#'
#' @param block.results a list containing two lists, one for finngen and one for ukbb, with the loaded RData in each entry
#'Note that here, the finngen ones are the ones changing each time, UKBB results stay the same.
#' @return list of lists containing the matching statistics for each pair gamma parameters
#' @export
#'
#' @examples
fullPairwiseComparisons <- function(block.results, len.run=11)
{
  reordered.lists.fg <- list()
  reordered.lists.uk <- list()
  stopifnot(length(block.results$finngen) == length(block.results$ukbb))
  adj.length = len.run-1
  for(i in 1:adj.length)
  {
    reordered.lists.fg[[i]] <- lapply(i:(adj.length+i), function(j) block.results$finngen[[(j%%len.run+1)]])
    names(reordered.lists.fg[[i]]) <- sapply(i:(adj.length+i), function(j) names(block.results$finngen)[(j%%len.run+1)])

    reordered.lists.uk[[i]] <- lapply(i:(adj.length+i), function(j) block.results$ukbb[[(j%%len.run+1)]])
    names(reordered.lists.uk[[i]]) <- sapply(i:(adj.length+i), function(j) names(block.results$ukbb)[(j%%len.run+1)])
  }
  #Need to also add the matching case to get the diagonal elements:
  reordered.lists.fg[[len.run]] <- block.results$finngen
  reordered.lists.uk[[len.run]] <- block.results$ukbb
  lapply(reordered.lists.fg, function(x) makeTableOfComparisons(block.results$ukbb, x,by_order = TRUE))
  #This "by_order" flag super important!

}


#' Function to take all the pairwise run data and put it into a table
#'
#' @param all.pairwise.settings a list of objects looking at the U and V comparison for pairiwise test
#' elements have U,V, Xhat norms, and more complete U and V data (distribution) [.e.g a list of objects outputted from makeTableOfComparisons]
#' @param block_results the results from both UKBB and finngen- the UKBB is used as the "unchanging" name reference list
#' @param changing.lists the lists of reordered gammas for a given sample, to match the "Finngen: reorderin,
#' since that one had the rotating gammas in fullPairwiseComparisons above (different gamma references)
#' @param pattern the regualr expression of gamma terms to extract (group 1)
#'
#' @return
#' @export
#'
#' @examples
extractAndTabularizeRun <- function(all.pairwise.settings, block_results, changing.lists, pattern = "adjusted_([0-9\\.]+)-shrinkage")
{
  growing.V.crosscohort <- NULL
  growing.U.crosscohort <- NULL
  Xhat_norms <- NULL
  for(i in 1:length(all.pairwise.settings))
  {
    growing.V.crosscohort <- rbind(growing.V.crosscohort, all.pairwise.settings[[i]]$V %>%
                                     mutate("UKBB.factors" = names(block_results$ukbb), "FG.factors"= names(changing.lists[[i]])) %>%
                                     mutate("Gamma.ukbb" = as.factor(seq(0,1,0.1))))
    growing.U.crosscohort <- rbind(growing.U.crosscohort, all.pairwise.settings[[i]]$U %>%
                                     mutate("UKBB.factors" = names(block_results$ukbb), "FG.factors"= names(changing.lists[[i]]))%>%
                                     mutate("Gamma.ukbb" = as.factor(seq(0,1,0.1))))
    Xhat_norms <- rbind(Xhat_norms, data.frame("dist" = unlist(all.pairwise.settings[[i]]$X_hat.norms),
                                               "Gamma.fg"=as.factor(stringr::str_match(names(changing.lists[[i]]), pattern)[,2]),  "Gamma.ukbb"=as.factor(seq(0,1,0.1))))
  }
  growing.V.crosscohort$Gamma.fg <- as.factor(stringr::str_match(growing.V.crosscohort$FG.factors, pattern)[,2])
  growing.U.crosscohort$Gamma.fg <- as.factor(stringr::str_match(growing.U.crosscohort$FG.factors, pattern)[,2])
  #Make all the relevant columns numeric:
  growing.V.crosscohort <- growing.V.crosscohort %>%
    mutate_at(c("correlation","PRC","GlobalKappa", "Kendall_correlation","Procrustes_pearson","Euclidian_dist"), as.numeric)
  growing.U.crosscohort <- growing.U.crosscohort %>%
    mutate_at(c("correlation","PRC","GlobalKappa", "Kendall_correlation","Procrustes_pearson","Euclidian_dist"), as.numeric)
  return(list("V" = growing.V.crosscohort %>% select(-Factors),
              "U"=growing.U.crosscohort %>% select(-Factors), "Xhat_norms"=Xhat_norms))
}


#' Function to take all the pairwise run data and put it into a table
#'
#' @param all.pairwise.settings a list of objects looking at the U and V comparison for pairiwise test
#' elements have U,V, Xhat norms, and more complete U and V data (distribution) [.e.g a list of objects outputted from makeTableOfComparisons]
#' @param block_results the results from both UKBB and finngen- the UKBB is used as the "unchanging" name reference list
#' @param changing.lists the lists of reordered gammas for a given sample, to match the "Finngen: reorderin,
#' since that one had the rotating gammas in fullPairwiseComparisons above (different gamma references)
#' @param pattern the regualr expression of gamma terms to extract (group 1)
#'
#' @return
#' @export
#'
#' @examples
extractAndTabularizeRunGeneric <- function(all.pairwise.settings)
{
  growing.V.crosscohort <- NULL
  growing.U.crosscohort <- NULL
  growing.paired.data <- NULL
  Xhat_norms <- NULL
  for(i in 1:length(all.pairwise.settings))
  {
    growing.V.crosscohort <- rbind(growing.V.crosscohort, all.pairwise.settings[[i]]$V)
    growing.U.crosscohort <- rbind(growing.U.crosscohort, all.pairwise.settings[[i]]$U)
    Xhat_norms <- rbind(Xhat_norms, data.frame("dist" = unlist(all.pairwise.settings[[i]]$X_hat.norms)))
    growing.paired.data <- rbind(growing.paired.data, all.pairwise.settings[[i]]$matrix_paired_results)
  }
  #Make all the relevant columns numeric:
  growing.V.crosscohort <- growing.V.crosscohort %>%
    mutate_at(colnames(growing.V.crosscohort)[4:17], as.numeric)
  growing.U.crosscohort <- growing.U.crosscohort %>%
    mutate_at(colnames(growing.V.crosscohort)[4:17], as.numeric)
  return(list("V" = growing.V.crosscohort,
              "U"=growing.U.crosscohort, "Xhat_norms"=Xhat_norms, "joint"=growing.paired.data))
}



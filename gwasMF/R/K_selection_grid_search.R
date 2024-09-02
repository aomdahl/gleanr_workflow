
#' Title
#'
#' @param opath
#' @param option
#' @param X
#' @param W
#' @param W_c
#' @param all_ids
#' @param names
#' @param step.limit - what is the maximum number of grid look steps you are willing to take?
#' @param init.range - how broad should the initial search range be? Default is 4
#' @param parallel - options are "none" or "multicore" or "snow"
#' @param ncpu - "how many CPUs to use if going parallel
#' @param ...
#'
#' @return the bic_dat object with the best BIC
#' @export
#'
#' @examples
gridSearchK <- function(opath, option,X,W,W_c,all_ids,names,step.limit = 10,init.range=4,...)
{
  #parallel settings
  parallel="no"
  if(option$ncores > 1){parallel = "multicore"; message("Running in parallel on ", option$ncores, " cores.")}
  #Specify a valid search range
  M = ncol(X)
  k.range <- 1:(M-1)

  #Special case- k.range <= 4:
  #Choose init params, 4:
  init.params <- list("K" = floor(quantile(k.range,seq(0, 1, 1/init.range)))[-1])

  #Edge case- make sure no bad initializations
  if(any(init.params$K <= 0))
  {
    d <- which(init.params$K ==0)
    init.params$K<- init.params$K[-d]
  }
  #Initial grid search
  #optimizeK <- function(K, opath, option, X_in, W_in, W_c, all_ids, names, reg.vect,...)
  init.results <- gridSearch(as.list(init.params$K), option$ncores, opath, option, X, W, W_c, all_ids, names,...)

  #init.results <- paramtest::grid_search(optimizeK, params = init.params, n.iter = 1 , boot = FALSE, bootParams = NULL, parallel = parallel, ncpus = option$ncores,
  #                                       opath = opath,
  #                                       option = option,
  #                                       X_in = X,
  #                                       W_in = W,
  #                                       W_c = W_c,
  #                                       all_ids = all_ids,
  #                                       names = names)


  sdv <- gatherSearchData(init.results,k.range)
  grid.search.record <- gridSearchRecord(init.results, init.params, NULL)
  if(sdv$next_params$terminate)
  {
    message("All BICs the same. Grid search is complete.")
    if(all(is.infinite(sdv$query_matrix$BIC)))
    {
      message("All entries are 0'd out. Program will end.")
    }
    return(sdv$min_result)
  }
  save(sdv,grid.search.record, file = paste0(opath, "_K_search.RData"))
  #If we aren't looking at immediatley adjacent K and we haven't exceeded our step limit
  while(length(sdv$query_matrix$K) < step.limit & sdv$k_diff >= 1)
  {
    #Keep looking
    next.grid <- gridSearch(as.list(sdv$next_params$K), option$ncores, opath, option, X, W, W_c, all_ids, names,...)
    #next.grid <- paramtest::grid_search(optimizeK, params = sdv$next_params, n.iter = 1 , boot = FALSE, bootParams = NULL,
    #                                    parallel = parallel, ncpus =option$ncores,
    #                                    opath = opath,
    #                                    option = option,
    #                                    X_in = X,
    #                                    W_in = W,
    #                                    W_c = W_c,
    #                                    all_ids = all_ids,
    #                                    names = names)
    #New paradigm- we need this to update off the full list
    grid.search.record <- gridSearchRecord(next.grid, sdv$next_params, grid.search.record) #use this last entry in the search
    sdv <- gatherSearchData(next.grid,k.range, curr_grid=sdv$query_matrix,curr_best = sdv$min_result ) #need previous dat here: UPDATE
    if(sdv$next_params$terminate)
    {
      message("All BICs the same. Grid search is complete.")
      if(all(is.infinite(sdv$query_matrix$BIC)))
      {
        message("All entries are 0'd out. Program will end.")
        return(sdv$min_result)
      }
    }

    save(sdv,grid.search.record, file = paste0(opath, "_K_search.RData"))
  #query.matrix <- storeGridResults(next.grid, curr_grid=query.matrix)
  #min.result <- getBestRun(next.grid,curr_best = min.result)
  #For the next iteration:
  #curr_best_K = query.matrix$K[which.min(query.matrix$BIC)]
  #next.params <- chooseNextParams(curr_best_K,query.matrix,k.range)
  #check we aren't looking at super close things
  #k.diff = checkK(c(curr_best_K,next.params$K))
  }
  #save(sdv, file="/scratch16/abattle4/ashton/DEBUG.RData")
  save(sdv,grid.search.record, file = paste0(opath, "_K_search.RData"))
  sdv$min_result
}

#' Manages results of grid search and nominates next parameters to look at
#'
#' @param gs_object - output of grid_search function
#' @param k_range - full range of possible K to consider searching
#' @param ...
#' May be curr_grid=sdv$query_matrix, the current matrix of searches, or
#' curr_best = sdv$query_matrix$min_result, the best bic_dat object so far
#' @return a list containing
#' @export
#'
#' @examples
#' #gatherSearchData(next.grid,k.range, curr_grid=sdv$query_matrix,curr_best = sdv$min_result )
gatherSearchData <- function(gs_object,k_range,curr_grid=NULL,curr_best=NULL)
{
  #Store results so we don't repeat tests
  query.matrix <- storeGridResults(gs_object,curr_grid = curr_grid)
  #Which one is best
  min.result <- getBestRun(gs_object,curr_best= curr_best)
  curr.best.K = query.matrix$K[which.min(query.matrix$BIC)]

  #What parameters should we consider next? Binary search approach
  next.params <- chooseNextParams(curr.best.K,query.matrix,k_range)
  #How close are the next Ks we are considering to existing Ks? Make sure we don't repeat anything.
  k.diff = checkK(next.params$K)
  return(list("query_matrix"=query.matrix, "min_result"=min.result,
              "curr_best_K"=curr.best.K,
              "next_params" = next.params, "k_diff"=k.diff))
}

#' Helper function to track the results of run from one to the nnext test
#'
#' @param gs_object return from grid search object
#' @param params which settings of K we tried at
#' @param record_obj the object this is being stored in
#'
#' @return
#' @export
#'
#' @examples
gridSearchRecord <- function(gs_object,prev.results, params, record_obj)
{
  if(is.null(record_obj))
  {
    record_obj <- list("test_dat"=list(), "curr_runs"=NULL)
  }
  bics <- sapply(gs_object$results, function(x) (x$min.dat$min_sum))
  #bics_global procedure
  #Get the scalar that is max across all for the variance comparison; generally should be M-1.
  bics_global <- getGlobalBICSum(gs_object$results,record_obj$test_dat) #need to include the previous information

  alphas <- sapply(gs_object$results, function(x) (x$min.dat$alpha))
  lambdas <- sapply(gs_object$results, function(x) (x$min.dat$lambda))
  K_end <- sapply(gs_object$results, function(x) ncol(x$optimal.v))
  K_end[is.null(K_end)] <- 0

  record_obj$curr_runs <- rbind(record_obj$curr_runs,
                                data.frame("Kinit"=params$K,
                                           "K_end"=K_end,
                                           "bic_local"=bics,
                                           "bic_global"=bics_global,
                                           "alpha"=alphas,
                                           "lambda"=lambdas))
  record_obj$test_dat <- c(record_obj$test_dat, gs_object) #concatenate all the new information together
  record_obj
}

#Only works for N =2 at the moment

#' Figure out which points to test next
#' Currently limited to BINARY SEARCH- only searches 2 points next, at the midpoint between points that have been previously tested (integers only)
#' @param best_K- the K that is currently best
#' @param curr_grid-a list with all the BICs and all the Ks
#' @param all_K- all possible Ks to consider
#'
#' @return
#' @export
#'#For debugging, try:
#load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/run_scripts/finngen_ukbb_benchmark/v2_expanded_run/bic_method_comparisons_covar_adjusted/GRID.DEBUG.RData")
#' @examples
#' #chooseNextParams(curr.best.K,query.matrix,k_range)
chooseNextParams <- function(best_K, curr_grid,all_K)
{
  end = FALSE
  if(all(curr_grid$BIC == curr_grid$BIC[which(curr_grid$K == best_K)]))
    {
      #All the BICs are the same, stop.
      end = TRUE
  }
  testbounds <- getNewTestBounds(best_K, curr_grid,all_K)
  above <- testbounds$upper; below <- testbounds$lower

    #A number in the range above and below:
    above.test <- floor(quantile(best_K:above,probs = c(0.5)))
    below.test <- floor(quantile(below:best_K,probs = c(0.5)))
    test.list <- c(above.test, below.test)

  #Ensure there are no duplicates in new list:
  if(any(duplicated(test.list)))
  {
    message("duplicates in list- not good.")
    test.list <- test.list[!duplicated(test.list)] #Drop duplicates
  }
  #Ensure we aren't repeating previous tests
  if(any(test.list %in% curr_grid$K))
  {
    message("Caught a repeat test, avoiding")
    repeats <- which(test.list %in% curr_grid$K)
    test.list <- test.list[-repeats]
  }

  list("K"=test.list, "terminate"=end)

}

#' Helper function for chooseNextParams to get the bounds of the space to test next.
#'
#' @param best_K
#' @param curr_grid
#' @param all_K
#'
#' @return upper (above) and lower (below)
#' @export
#'
#' @examples
#' #getNewTestBounds(best_K, curr_grid,all_K)
getNewTestBounds <- function(best_K, curr_grid,all_K)
{
  diffs <- (curr_grid$K - best_K)
  if(best_K == max(curr_grid$K))
  {
    above <- max(all_K)
    #Special case- we are on the edge of what we can try
    if(best_K == max(all_K))  {    above <- max(all_K) - 1    }
    below <-  curr_grid$K[which(diffs == max(diffs[diffs < 0]))]
  } else if(best_K == min(curr_grid$K))
  {
    above <- curr_grid$K[which(diffs == min(diffs[diffs > 0]))]
    below <- min(all_K)
    #Special case- we are on the edge of what we can try
    if(best_K == min(all_K))  {    below <- min(all_K) + 1    }
  }else
  {
    #closest above (omitting the 0 case)
    above <- curr_grid$K[which(diffs == min(diffs[diffs > 0]))]
    below <-  curr_grid$K[which(diffs == max(diffs[diffs < 0]))]
  }
  list("upper"=above, "lower"=below)
}

#' Update search grid with newest search data.
#'
#' @param gs_object list of bic_autofit objects
#' @param curr_grid current data containing a list of Kinits and global BIC scores
#'
#' @return an updated grid with the latest looks.
#' @export
#'
#' @examples
storeGridResults <- function(gs_object, curr_grid=NULL)
{
  tested_k <- gs_object$tests[,2]
  bic_scores <- sapply(gs_object$results, function(x) x$min.dat$min_sum)
  bic_scores <- getGlobalBICSum(gs_object$results)
  ret <- list("K"=tested_k, "BIC"=bic_scores)
  if(!is.null(curr_grid))
  {
    ret <- curr_grid
    ret$K <- c(ret$K,tested_k)
    ret$BIC <- c(ret$BIC,bic_scores)
  }
  ret
}
#
getBestRun <- function(gs_object, curr_best = NULL)
{
  #$$todo
  # Instead of just looking at the best score, we need to rescale all to be on the same variance scale.
  #Do this by getting the fit.scalar for the object with the largest K across the set under consideration
  #then get new BIC scores for each with fit.term/global.fit.scalar + df.term + addends
  #Note that for sklearn, the addends will be problematic since they are subject to the sse of that calculation.
  #They also contain ebic
  #consider just dropping this
  #From all of these rescaled BIC terms, pick the one that minimizes
  #make sure to test in zou case:
  #we pick a terrible zou because it is scaled to look good even though its terrible.
  #list("bic.list" = BIC,
  #     "fit.term" = deviance(fit),
  #     "df.term"=  log(n)*k,
  #     "fit.scaler"=1,
  #     "addends" =0)

  #min_i <- which.min(sapply(gs_object$results, function(x) x$min.dat$min_sum))
  #Its possible that there are matching minimums.
  #In this case, choose the one that has the smaller Kinit score
  find_min_of <- getGlobalBICSum(gs_object$results)
  min_i = which(find_min_of == min(find_min_of, na.rm = TRUE))
  if(length(min_i) > 1)
  {
    min.k.i <- which.min(gs_object$tests$K[min_i])
    min_i= min_i[min.k.i]
  }
  #min_i <- which.min(getGlobalBICSum(gs_object$results))
  new.best <- gs_object$results[[min_i]]
  if(!is.null(curr_best))
  {
    if(curr_best$min.dat$min_sum < new.best$min.dat$min_sum) #if the current best is better than the new best
    {
      new.best <- curr_best #don't change the new best
      message("no change!")
    }
  }
  new.best
}

  checkK <- function(ks)
  {
    #Have to check, if we are testing the same K twice, need to stop
    if(any(duplicated(ks)))
    {
      #return(0)
      ks <- ks[!duplicated(ks)] #Drop duplicates
    }
    if(length(ks) == 0)
    {
      return(0)
    }
    if(length(ks) == 1) #We only have 1 to test
    {
      return(1)
    }
    #However, in the grid any item - itself is 0.
    tests <- expand.grid(ks,ks)
    nonzero.diffs <- abs(tests[,1] - tests[,2])
    z <- which(nonzero.diffs == 0)

    min(nonzero.diffs[-z])

  }
# Load the rgenoud package
# Define a wrapper function for genoud to minimize bic.sum by modifying K
optimizeK <- function(K, opath, option, X_in, W_in, W_c, all_ids, names,...) {
  option$K <- K
  bic.dat <- getBICWorkhorse(opath, option, X_in, W_in, W_c, all_ids, names, ...)
  message("K: ",K, ", BIC: ", bic.dat$min.dat$min_sum)
  return(bic.dat)
}


gridSearch <- function(iter_params, ncpus, opath, option, X_in, W_in, W_c, all_ids, names,...)
{
  ret.list <- list("tests"=data.frame("iter"=1:length(iter_params), "K"=unlist(iter_params)))
  if(ncpus == 1 | length(iter_params) == 1)
  {
    message("Not running in parallel this time...")
    ret.list$results <- lapply(iter_params, function(kinit) optimizeK(kinit, opath, option, X_in, W_in, W_c, all_ids, names,...))
  }else
  {
    message("Now running in parallel across ", ncpus, " cores.")
    ret.list$results <- parallel::mclapply(iter_params, function(kinit)
    {
      optimizeK(kinit, opath, option, X_in, W_in, W_c, all_ids, names,...)
    }, mc.cores = ncpus)
  }
  ret.list
}

getGlobalBICs <- function(fit_vect, source)
{
  #
  max.coef.index <- which.max(sapply(fit_vect, function(x) (x$min.dat[[source]]$n.coef)))
  scalar_global <- fit_vect[[max.coef.index]]$min.dat[[source]]$fit.scaler #
  #Get the scalar that is max across all for the variance comparison; generally should be M-1.
 sapply(fit_vect, function(x) x$min.dat[[source]]$fit.term/scalar_global + x$min.dat[[source]]$df.term + x$min.dat[[source]]$addends)
}


getGlobalBICSum <- function(fit_vect)
{
  getGlobalBICs(fit_vect, "bic_a_dat") + getGlobalBICs(fit_vect, "bic_l_dat")
}

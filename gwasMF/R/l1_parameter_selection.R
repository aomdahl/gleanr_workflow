#Parameter selection helper functions
singlePointGridExpansion <- function(sparsity.params, optimal.sparsity.param, n.points)
{
  sorted.sparsity.params <- sort(sparsity.params, index.return = TRUE)
  ordered.list <- sorted.sparsity.params$x
  optimal.index = which.min(abs(optimal.sparsity.param-sparsity.params))
  sorted.index <- which(sorted.sparsity.params$ix == optimal.index)
  if(min(ordered.list) == optimal.sparsity.param)
  {
    above <- ordered.list[sorted.index + 1]
    below <- optimal.sparsity.param - (above - optimal.sparsity.param)
    if(below < 0)
    {
      below = 0
      lower = optimal.sparsity.param * (1/(n.points - 1))
      #seq(0,optimal.sparsity.param, length.out=n.points-1 )
    }
  } else if(max(ordered.list) == optimal.sparsity.param) #its the largest parameter tested
  {
    below <- ordered.list[sorted.index - 1]
    above <- optimal.sparsity.param + (optimal.sparsity.param-below)
  }else {
    #Its bounded- our estimates should be between the one immediately above and below
    above <- ordered.list[sorted.index + 1]
    below <- ordered.list[sorted.index - 1]
  }
  if(is.na(above) | is.na(below))
  {
    print("proposed new paramters are not possible")
    quit()
  }
  if(above == below)
  {
    message("Converged on single solution")
    return(below)
  }
  new.list <- seq((below),(above),length.out=n.points+2)
  new.list <- new.list[-c(1, length(new.list))]
  if(any(new.list <= 0))
  {
    print("removing negatives...")
    rep.list <- new.list[which(new.list > 0)]
    if(any(new.list == 0))
    {
      rep.list <- c(min(rep.list)/2,rep.list)
    }
    new.list <- rep.list
  }
  new.list
}

#we don't want to be redundant, so omit the ones that are repeated...
#alphas: the optimal parameters to consisder
#all.a: the full list of parameters from which they are drawn
#n.points: how many new points to pick from.
ProposeSparsityParamsFromGrid <- function(alphas,all.a, n.points)
{
  #N-points will include the original ones...
  if(length(alphas) == 1)
  {
    new.list <- singlePointGridExpansion(all.a[order(all.a)], alphas, n.points)
  }else
  {
    range.a <- range(alphas)
    new.list <- (seq(range.a[1],range.a[2],length = n.points))
  }
  if(any(is.na(new.list)))
  {
    message("detected possible ERROR in proposed parameters")
    print(new.list)
    message("Dropping NAs")
    new.list <- new.list[!is.na(new.list)]
  }
  if(any(new.list < 0))
  {
    message('Detected negatives in proposed parameters; dropping these')
    new.list <- new.list[new.list > 0]
  }
  if(length(new.list) == 0)
  {
    message("Possible error, no terms detected")
    message("Retruning original input.")
    return(alphas)
  }
return(new.list)

}


ProposeNewSparsityParamsFromScore <- function(bic.list,sparsity.params, n.points = 7, no.score = FALSE)
{
  #we want to keep it griddy. So a bit further downstream compare this to the BIC, I think this actually the one that we want, but okay.
  if(no.score)
  {
    #then bic.list is the optimal one; generate fake scores
    fake.scores <- rep(100,length(sparsity.params))
    fake.scores[which(sparsity.params == bic.list)] <- -1
    bic.list <- fake.scores
  }
  optimal.sparsity.param <- sparsity.params[which.min(bic.list)]
  message("optimal sparsity param is ", optimal.sparsity.param)
  ordered.list <- sparsity.params[order(sparsity.params)] #order the list
  #New paradigm: always look above and below,
  #If its the smallest paramter tested
  singlePointGridExpansion(sparsity.params[order(sparsity.params)], optimal.sparsity.param, n.points)
}


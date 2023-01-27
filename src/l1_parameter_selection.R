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

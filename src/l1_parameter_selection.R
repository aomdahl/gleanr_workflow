#Parameter selection helper functions
singlePointGridExpansion <- function(ordered.list, optimal.sparsity.param, n.points)
{
  if(min(ordered.list) == optimal.sparsity.param)
  {
    above <- ordered.list[which(ordered.list == optimal.sparsity.param) + 1]
    below <- optimal.sparsity.param - (above - optimal.sparsity.param)
    if(below < 0)
    {
      below = 0
      lower = optimal.sparsity.param * (1/(n.points - 1))
      #seq(0,optimal.sparsity.param, length.out=n.points-1 )
    }
  } else if(max(ordered.list) == optimal.sparsity.param) #its the largest parameter tested
  {
    below <- ordered.list[which(ordered.list == optimal.sparsity.param) - 1]
    above <- optimal.sparsity.param + (optimal.sparsity.param-below)
  }else {
    #Its bounded- our estimates should be between the one immediately above and below
    above <- ordered.list[which(ordered.list == optimal.sparsity.param) + 1]
    below <- ordered.list[which(ordered.list == optimal.sparsity.param) - 1]
  }
  l=as.integer((n.points)/2) + 2
  upper <- seq(optimal.sparsity.param, above, length.out = l)
  lower <- seq(below, optimal.sparsity.param, length.out = l)
  new.list <- c(lower[2:(l-1)], optimal.sparsity.param, upper[2:(l-1)])
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

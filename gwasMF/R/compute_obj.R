##############################################################################################################
## Compute objective: ||(X - LF') .* W||_F^2 + alpha1*|L|_1 + lambda1*|F|_1
##############################################################################################################

compute_obj <- function(X, W, L, FactorM, option, decomp = FALSE){
	
	if(is.null(FactorM) | is.null(L))
	{
		message("Unable to calculate objective for this run, because no valid matrix generated")
		return(NA)
	}
  if(is.null(option$alpha1) | is.null(option$lambda1))
  {
    message("Unable to calculate objective, no valid sparsity parameters given.")
    return(NA)
  }
  M = ncol(X)
  N = nrow(X)
  residuals = (W * X - W * (L %*% t(FactorM)));
	#Residual_penalty = sum(sum(residual^2));
#	Residual_penalty = sum(residual^2) / (2 * nrow(L));
#After some quick evaluation, we either scale everything or don't scale at all. For simplicity we scale nothing.
	
	# the l2 penalty is added in the fitting step to avoid singularity, and is accounted for here
	L_penalty = option[['alpha1']]* sum(abs(L)) + 1e-20 * sqrt(sum(L^2));
	FactorM_penalty = option[['lambda1']] * sum(abs(FactorM)) + 1e-20 * sqrt(sum(FactorM ^ 2));
  if(option$fixed_ubiq)
  {
    FactorM_penalty = option[['lambda1']] * sum(abs(FactorM[,-1])) + 1e-20 * sqrt(sum(FactorM[,-1] ^ 2));
   
  }

  #derived version
  Residual_penalty = (M+N)/(2*M*N)*sum(residuals^2)
  obj = Residual_penalty + FactorM_penalty + L_penalty;   
  #AVerage version
  avg_v  = (1/M) * ((1/2*N) * sum(residuals)^2 + FactorM_penalty)
  avg_u = (1/N)  * ((1/2*M) * sum(residuals)^2 + L_penalty)                                                                      
  #print(paste0('obj decomposition: ', (Residual_penalty), '; ', (L_penalty), '; ', (FactorM_penalty)));
  #print(paste0("Average terms: V: ", avg_v, "  U: ", avg_u))
  #reportObjectivePerFactor(X,W,L,FactorM, option)
  if(decomp)
  {
    return(list("Fit.p" = Residual_penalty, "V.p"=FactorM_penalty, "U.p"= L_penalty, "Tot" = obj))
  }
	return(obj)
}
reportObjectivePerFactor <- function(X,W,L,FactorM, option)
{
  init = 1
  if(option$fixed_ubiq)
  {
    init = 2
  }
  #Reorder factors by PVE
  pve.rank=order(PercentVarEx(as.matrix(X)*as.matrix(W), v = V), decreasing = TRUE)
  curr.indices <- c()
  for(i in init:ncol(FactorM))
  {
    index = which(pve.rank == i)
    curr.indices <- c(curr.indices, index)
    print(curr.indices)
    #does adding this in increase our objective
    #What do we do with the objective and ubiq 1?
    print(compute_obj(X, W, L[,curr.indices], FactorM[,curr.indices], option))
  }
  
}
#option <- list()
#option$alpha1 <- f.qd$autofit_alpha[1]
#option$lambda1 <- f.qd$autofit_lambda[1]
#option$fixed_ubiq <- FALSE


#  Residual_penalty = (M+N)/(2*M*N)*sum(residuals^2)
#residuals = W * (X - L %*% t(FactorM));

calcMSE <- function(X,W,Fa, L)
{
  sum(X*W- (L%*% t(Fa))*W)^2
}

#Beta test function- can we pick out the best factors by how they contribute (or don't) to the objective?
#This will balance in sparsity better than PVE does.
#DropFactorsByObjective(X,W,reg.run$U,reg.run$V, minK=maxK, option, maxK = maxK)
DropFactorsByObjective <- function(X,W,U,V, minK, option, maxK = NULL,drop.min.change = TRUE)
{
	if(is.null(maxK) | maxK == 0)
	{
		maxK = ncol(X) - 1
	}
  #1) End if only 1 column left doing fixed ubiq version
  if(ncol(as.matrix(V)) == 1 && option$fixed_ubiq)
  {
	  return(list("U"=U, "V" = V, "K" = ncol(V)))
  }
  #2) End if there are no more columns left to remove
  if(ncol(as.matrix(V)) < 1)
  {
	  message("All factors removed; none contribute to objective")
	  return(list("U"=U, "V" = V, "K" = NA))
  }
  remaining.set <- 1:ncol(V)
  #Repeat the process until conditions:
  #we reach our minimum k size threshold
  if(length(remaining.set) <=  minK)
  {
    return(list("U"=U, "V" = V, "K" = ncol(V)))
  }
  # function(X, W, L, FactorM, option)
  init.obj <- compute_obj(X,W, U, V, option)
  drop.set <- remaining.set
  if(option$fixed_ubiq)
  {
    drop.set <- 2:ncol(U)
  }
  min.obj <- init.obj
  min.index <- -1
  #Try dropping each one
  #if removing one reduces the objective, drop that one that results in biggest reduction.
  print("new recursive level...")
  print("init.obj.")
  print(init.obj)
  #A greedy approach to remove factors- drop the one that minimizes the objective overall
  for(i in drop.set) #we don't drop ubiq one
  {
    new.obj <- compute_obj(X,W, U[,-i], V[,-i], option)
     if(new.obj < min.obj)
    {
      min.index <- i
      min.obj <- new.obj
    }
  }
  
  #Special case- continue to remove until reached a number.∂f
  if(drop.min.change & min.index == -1 & ncol(V) < maxK)
  {
    message("Continuing to remove factors, even though removing more hurts the objective")
    #Iterate through and get the one with the lowest drop
    min.increase = Inf
    for(i in drop.set) #we don't drop ubiq one
    {
      new.obj <- compute_obj(X,W, U[,-i], V[,-i], option)
      diff = new.obj - init.obj
      if(diff < min.increase)
      {
        min.index <- i
        min.increase <- diff
      }
    }
    return(DropFactorsByObjective(X,W,U[,-min.index],V[,-min.index], minK, option, maxK = maxK, drop.min.change = TRUE))
  }
  
  
  #Removing any objective increases (or stays the same) rather than decreases.
  if((min.obj == init.obj | min.index == -1) & (ncol(U) > maxK))
  {
    #Now we need to call the funciton again, except this time specify the condition of choosing the one which increases the objective the smallest amount.
    return(DropFactorsByObjective(X,W,U[,-min.index],V[,-min.index], minK, option, maxK = maxK, drop.min.change = TRUE))
  }
  #3) End if dropping any more factors increases the objective and we have reached the max number of K 
  else if((min.obj == init.obj | min.index == -1) & (ncol(U) <= maxK))
  {
    return(list("U"=U, "V" = V, "K" = ncol(V)))
  }else
  {
    return(DropFactorsByObjective(X,W,U[,-min.index],V[,-min.index], minK, option))
  }
}

PruneFactorsByObjective <- function(X,W,U,V, minK, option)
{
  keep.dat <- DropFactorsByObjective(X,W,U,V, minK, option)
  return(list("U"))
}

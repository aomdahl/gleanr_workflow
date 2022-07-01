##Helper regression functions for organization.
  suppressWarnings(library(penalized))
  library(foreach)
  library(doParallel)
 glmnetLASSO <- function(dat_i,xp, colcount, lambdas, penalties)
{
    out <- tryCatch(
    expr = { penalties <- c(0, rep(1, (colcount - 1)))
      fit <- glmnet(x = dat_i[,paste0('F', seq(1,colcount))], y = xp, alpha = 1, lambda = lambdas,
                    intercept = FALSE, penalty.factor = penalties)
      coef(fit, 'all')[,1][-1];
        
    },
    error = function(e){
	print(e)
        print("Vals")
        print(head(dat_i[,paste0('F', seq(1,colcount))]))
        print("Outcomes")
        print(xp)
	readline()
        rep(0,colcount)

    },
     warning = function(w){
            message('Caught an warning!')
            print(w)})
    return(out)
}
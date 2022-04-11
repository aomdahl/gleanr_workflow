penalizedDEBUG <- function (response, penalized, unpenalized, lambda1 = 0, lambda2 = 0, 
          positive = FALSE, data, fusedl = FALSE, model = c("cox", 
                                                            "logistic", "linear", "poisson"), startbeta, startgamma, 
          steps = 1, epsilon = 1e-10, maxiter, standardize = FALSE, 
          trace = FALSE) 
{
  message("here")
  if (missing(maxiter)) 
    maxiter <- if (lambda1 == 0 && lambda2 == 0 && !positive) 
      25
  else Inf
  if (steps == "Park" || steps == "park") {
    steps <- 1
    park <- TRUE
  }
  else park <- FALSE
  prep <- penalized:::.checkinput(match.call(), parent.frame())
  if (ncol(prep$X) >= nrow(prep$X) && all(lambda1 == 0) && 
      all(lambda2 == 0) && !any(prep$positive)) 
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.", 
         call. = FALSE)
  fit <- penalized:::.modelswitch(prep$model, prep$response, prep$offset, 
                      prep$strata)$fit
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  nr <- nrow(prep$X)
  fusedl <- prep$fusedl
  if (length(lambda1) == pp && (!all(lambda1 == 0))) {
    wl1 <- c(numeric(pu), lambda1)
    lambda1 <- 1
  }
  else {
    wl1 <- 1
  }
  if (length(lambda2) == pp) 
    lambda2 <- c(numeric(pu), lambda2)
  print(park)
  print(steps)
  if (park || steps > 1 && fusedl == FALSE) {
    if (pu > 0) 
      lp <- drop(prep$X[, 1:pu, drop = FALSE] %*% prep$nullgamma)
    else lp <- numeric(n)
    chck <- (wl1 > 0) & c(rep(FALSE, pu), rep(TRUE, pp))
    gradient <- drop(crossprod(prep$X[, chck, drop = FALSE], 
                               fit(lp)$residuals))
    print(gradient)
    if (length(wl1) > 1) {
      rel <- gradient/(wl1[chck] * prep$baselambda1[chck])
    }
    else {
      rel <- gradient/(wl1 * prep$baselambda1[chck])
    }
    from <- max(ifelse(prep$positive[chck], rel, abs(rel)))
    message("Max is:")
    message(from)
    #readline()
    if (from < lambda1) {
      warning("Chosen lambda1 greater than maximal lambda1: \"steps\" argument ignored")
      steps <- 1
      park <- FALSE
      from <- lambda1
    }
  }
  else {
    message("didn't go...")
    from <- lambda1
  }
  lambda1s <- seq(from, lambda1, length.out = steps)
  beta <- prep$beta
  louts <- if (park) 
    4 * pp
  else length(lambda1s)
  outs <- vector("list", louts)
  rellambda1 <- lambda1s[1]
  ready <- FALSE
  i <- 0
  while (!ready) {
    ready <- (rellambda1 == lambda1)
    i <- i + 1
    if (!fusedl) {
      if (rellambda1 != 0 || any(prep$positive)) {
        if (all(lambda2 == 0)) {
          out <- penalized:::.steplasso(beta = beta, lambda = rellambda1 * 
                              wl1 * prep$baselambda1, lambda2 = 0, positive = prep$positive, 
                            X = prep$X, fit = fit, trace = trace, epsilon = epsilon, 
                            maxiter = maxiter)
        }
        else {
          out <-  penalized:::.lasso(beta = beta, lambda = rellambda1 * 
                          wl1 * prep$baselambda1, lambda2 = lambda2 * 
                          prep$baselambda2, positive = prep$positive, 
                        X = prep$X, fit = fit, trace = trace, epsilon = epsilon, 
                        maxiter = maxiter)
        }
      }
      else {
        if (pp > n) {
          P <-  penalized:::.makeP(prep$X, lambda2 * prep$baselambda2)
          gams <-  penalized:::.solve(crossprod(t(P)), P %*% beta)
          PX <- P %*% t(prep$X)
          Pl <- P * matrix(sqrt(lambda2 * prep$baselambda2), 
                           nrow(P), ncol(P), byrow = TRUE)
          PlP <- crossprod(t(Pl))
          out <-  penalized:::.ridge(beta = gams, Lambda = PlP, X = t(PX), 
                        fit = fit, trace = trace, epsilon = epsilon, 
                        maxiter = maxiter)
          out$beta <- drop(crossprod(P, out$beta))
        }
        else {
          out <-  penalized:::.ridge(beta = beta, Lambda = lambda2 * 
                          prep$baselambda2, X = prep$X, fit = fit, 
                        trace = trace, epsilon = epsilon, maxiter = maxiter)
        }
      }
    }
    if (fusedl) {
      out <-  penalized:::.flasso(beta = beta, lambda1 = rellambda1 * 
                       wl1 * prep$baselambda1, lambda2 = lambda2 * prep$baselambda2, 
                     chr = prep$chr, positive = prep$positive, X = prep$X, 
                     fit = fit, trace = trace, epsilon = epsilon, 
                     maxiter = maxiter)
    }
    if (trace) 
      cat("\n")
    beta <- out$beta
    if (!ready) {
      if (!fusedl) {
        if (park) {
          newpark <-  penalized:::.park(beta = beta, lambda = rellambda1 * 
                             wl1 * prep$baselambda1, lambda2 = 0, positive = prep$positive, 
                           X = prep$X, fit = out$fit)
          rellambda1 <- rellambda1 * (1 - newpark$hh)
          if (rellambda1 < lambda1 || rellambda1 == Inf) {
            rellambda1 <- lambda1
            beta <- out$beta
          }
          else {
            beta <- newpark$beta
          }
          lambda1s <- c(lambda1s, rellambda1)
        }
        else {
          rellambda1 <- lambda1s[i + 1]
          beta <- out$beta
        }
      }
      else {
        rellambda1 <- lambda1s[i + 1]
        beta <- out$beta
      }
    }
    outs[[i]] <- out
  }
  if (length(lambda2) > 1) 
    lambda2 <- lambda2[pu + 1:pp]
  outs <- sapply(1:i, function(nr) {
    thislambda1 <- lambda1s[[nr]] * ifelse(length(wl1) > 
                                             1, wl1[pu + 1:pp], wl1)
    penalized:::.makepenfit(outs[[nr]], pu, fusedl = fusedl, prep$model, 
                thislambda1, lambda2, prep$orthogonalizer, prep$weights, 
                prep$formula, rownames(prep$X))
  })
  if (length(outs) == 1) 
    outs <- outs[[1]]
  outs
}


#Range recommender
recommendRange <- function (response, penalized, unpenalized, lambda1 = 100, lambda2 = 0, 
                         positive = FALSE, data, fusedl = FALSE, startbeta, startgamma, 
                         steps = 1, epsilon = 1e-10, maxiter) 
{
  trace <- FALSE
  standardize <- FALSE
  park <- FALSE
  steps = 100
  if (missing(maxiter)) 
    maxiter <- if (lambda1 == 0 && lambda2 == 0 && !positive) 
      25
  else Inf
  prep <- penalized:::.checkinput(match.call(), parent.frame())
  if (ncol(prep$X) >= nrow(prep$X) && all(lambda1 == 0) && 
      all(lambda2 == 0) && !any(prep$positive)) 
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.", 
         call. = FALSE)
  fit <- penalized:::.modelswitch(prep$model, prep$response, prep$offset, 
                                  prep$strata)$fit
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  nr <- nrow(prep$X)
  fusedl <- prep$fusedl
  if (length(lambda1) == pp && (!all(lambda1 == 0))) {
    wl1 <- c(numeric(pu), lambda1)
    lambda1 <- 1
  }
  else {
    wl1 <- 1
  }

  if (park || steps > 1 && fusedl == FALSE) {
    if (pu > 0) 
      lp <- drop(prep$X[, 1:pu, drop = FALSE] %*% prep$nullgamma)
    else lp <- numeric(n)
    chck <- (wl1 > 0) & c(rep(FALSE, pu), rep(TRUE, pp))
    gradient <- drop(crossprod(prep$X[, chck, drop = FALSE], 
                               fit(lp)$residuals))
    if (length(wl1) > 1) {
      rel <- gradient/(wl1[chck] * prep$baselambda1[chck])
    }
    else {
      rel <- gradient/(wl1 * prep$baselambda1[chck])
    }
    from <- max(ifelse(prep$positive[chck], rel, abs(rel)))
  }
   return(from)

}

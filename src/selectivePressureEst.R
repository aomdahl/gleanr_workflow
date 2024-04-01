########################################################################################################
## April 1, 2024
## Ashton Omdahl
## Assessing GLEANER
## Functions based on the work of SBayesR (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7896067/), also described in the supplement of 
## https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-020-0600-y/MediaObjects/41588_2020_600_MOESM1_ESM.pdf
## to estimate the relationship between MAF and effect size (something like selective pressur3) 
########################################################################################################

########################################### Simple base case ########################################### 
## (For testing) Simple case, 1 parameter- all non-zero effects, no sigma_g is assumed to be fixed
## do NOT recommend using this, mostly for testing

#' negative Logliklihood of normally distributed data, accounting for a selective effect (S) with fixed sigma
#'
#' @param s reflects relationship between MAF and effect size
#' @param sigma_g assumed fixed here, first temr omitted
#' @param effects GWAS or latent GWAS snp effects (non-zero effects only)
#' @param mafs vector of minor allele frequencies
#'
#' @return the negative log liklihood, for use with optim
#' @export
neglogLikFixedSigma <- function(s,sigma_g, effects, mafs)
{
  n <- length(effects)
  mafs.prod <- 2*mafs*(1-mafs)
  #-(-(n/2)* log(sigma_g) - (s*n)*sum(log(mafs.prod)) - sum((1/((mafs.prod)^s * sigma_g*2)) * effects^2))
  -(-(s/2)*sum(log(mafs.prod)) - sum((1/((mafs.prod)^s * sigma_g)) * effects^2))
}

#' Learn S using optim
#'
#' @param nonzero.betas vector of beta effect sizes, (non-zero effects only)
#' @param snp.mafs MAfs of those snps
#' @param sigma.g.fixed desired fixed sigma_g value
#'
#' @return optim output based on negative log liklihood
#' @export
estimateS_noSigmaG <- function(nonzero.betas, snp.mafs, sigma.g.fixed=1)
{
  optim(-0.001, fn = neglogLikFixedSigma,effects = nonzero.betas, mafs = snp.mafs, sigma.g.fixed = 1, method = "Brent", lower = -100, upper = 0 )
}
 
################################################# Slab-only case (0-effects omitted) ########################################### 
## Simple case- try to learn S (alpha in LDAK) and sigma_g^2, assuming the 0-effects have been dropped.

#' Negative log-liklihood of normally distributed data, with variance accounting for MAF, selective pressure, and study-wide variance
#'
#' @param par includes (s:a float reflecting relationship between MAF and effect size, sigma_g:scaling variance term)
#' @param effects GWAS or latent GWAS snp effects (non-zero effects only)
#' @param mafs cooresponding vector of minor allele frequencies
#'
#' @return the negative log liklihood, for use with optim
#' @export
neglogLik <- function(par, effects, mafs)
{
  s <- par[1]
  sigma_g <- par[2]
  n <- length(effects)
  mafs.prod <- 2*mafs*(1-mafs)
  -(-(n/2)* log(sigma_g) - (s/2)*sum(log(mafs.prod)) - sum((1/((mafs.prod)^s * sigma_g*2)) * effects^2))
}

#' Learn S and sigma_g^2 sing optim
#'
#' @param nonzero.betas vector of beta effect sizes, (non-zero effects only)
#' @param snp.mafs MAfs of those snps
#' @param init.s starting point for s search
#' @param init.sigma starting point for sigma search
#' 
#' @return optim output based on negative log liklihood
#' @export
estimateS_SigmaG <- function(nonzero.betas, snp.mafs, init.s = -0.001, init.sigma=0.01)
{
  optim(c(init.s, init.sigma), fn = neglogLik,effects = nonzero.betas, mafs = snp.mafs)
}

#################################################### Spike-and-slab case with unknown sparsity ########################################### 
## Harder case- learn S, sigma_g^2, and sparsity parameter pi

#' Negative log-liklihood of normally distributed data, with variance accounting for MAF, selective pressure, study-wide variance, and the proportion of 0-terms (pi)
#'
#' @param par includes (s:a float reflecting relationship between MAF and effect size, sigma_g:scaling variance term, pi:proportion of 0s)
#' @param effects GWAS or latent GWAS snp effects (non-zero effects only)
#' @param mafs cooresponding vector of minor allele frequencies
#'
#' @return the negative log liklihood, for use with optim
#' @export
neglogLikPi <- function(par, effects, mafs)
{
  s <- par[1]
  sigma_g <- par[2]
  pi <- par[3]
  n <- length(effects)
  mafs.prod <- 2*mafs*(1-mafs)
  -(n * log(pi) -(n/2)* log(sigma_g) - (s/2)*sum(log(mafs.prod)) - sum((1/((mafs.prod)^s * sigma_g*2)) * effects^2) + n*log(1-pi) + sum(ifelse(effects == 0, 1, 0)))
}

#' Learn S, sigma_g^2 and pi using optim
#'
#' @param nonzero.betas vector of beta effect sizes, (non-zero effects only)
#' @param snp.mafs MAfs of those snps
#' @param init.s starting point for s search
#' @param init.sigma starting point for sigma search
#' @param init.pi starting point for pi search (proportion of zero-sized effects.)
#' 
#' @return optim output based on negative log liklihood
#' @export
estimateS_SigmaG_pi <- function(all.betas, snp.mafs, init.s = -0.001, init.sigma=0.01, init.pi = 0.1)
{
  optim(c(init.s,init.sigma, init.pi), fn = neglogLikPi,effects = all.betas, mafs = snp.mafs)
}
#################################################### functions for unit testing ########################################### 
##Unit testing

## Calculate the variance term for a normal distribution based on f (MAF), s (selective pressure effect), and sigma_g (a global variance term.)
getVar <- function(f,s,sigma_g)
{
  ((2 * f * (1-f))^s) * sigma_g
}


testPiVersion <- function(seed)
{
  set.seed(seed)
  mafs <- runif(10000, min=0.001,max = 0.5)
  sigma_g = 0.04
  alpha = -0.9
  pi = 0.3
  all.vars <- sapply(mafs, function(f) getVar(f,alpha, sigma_g))
  print(summary(all.vars))
  sim.dat <- sapply(mafs, function(f)  rbinom(1, 1, pi) * rnorm(1, mean = 0, sd = sqrt(getVar(f,alpha, sigma_g))))
  print(hist(sim.dat, breaks = 30))
 optim(c(-0.001, 0.010, 0.01), fn = neglogLikPi,effects = sim.dat, mafs = mafs)
}


testSOnlyVersion <- function(seed)
{
  set.seed(seed)
  mafs <- runif(10000, min=0.001,max = 0.5)
  sigma_g = 0.04
  alpha = -0.9
  
  all.vars <- sapply(mafs, function(f) getVar(f,alpha, sigma_g))
  print(summary(all.vars))
  sim.dat <- sapply(mafs, function(f) rnorm(1, mean = 0, sd = sqrt(getVar(f,alpha, sigma_g))))
  hist(sim.dat, breaks = 30)             
  
  #PLot the negative log liklihood just to see......
  s.options <- seq(-0.001,-10,by=-0.005)
  find.val <- sapply( seq(-0.001,-10,by=-0.005), function(x) neglogLik(x,sigma_g,sim.dat, mafs))
  
  #This doesn't look right.
  print(plot(s.options, log(find.val)))
  
 optim(-0.001, fn = neglogLikFixedSigma,effects = sim.dat, mafs = mafs, sigma_g = sigma_g, method = "Brent", lower = -100, upper = 0 )
  
}

testS_SigmaVersion <- function(seed)
{
  set.seed(seed)
  mafs <- runif(10000, min=0.001,max = 0.5)
  sigma_g = 0.04
  alpha = -0.9
  
  all.vars <- sapply(mafs, function(f) getVar(f,alpha, sigma_g))
  print(summary(all.vars))
  print(summary(all.vars))
  sim.dat <- sapply(mafs, function(f) rnorm(1, mean = 0, sd = sqrt(getVar(f,alpha, sigma_g))))
  print(hist(sim.dat, breaks = 30))
  s.options <- seq(-0.001,-10,by=-0.05)
  sigma.g.options <- seq(0.001, 2, by = 0.01)
  optim(c(-0.001, 0.01), fn = neglogLik,effects = sim.dat, mafs = mafs)
}
#PLot the negative log liklihood just to see......


               


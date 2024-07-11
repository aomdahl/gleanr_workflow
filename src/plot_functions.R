#R plots data viz.
pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr, argparse, cowplot, ashr)

orderFactors <- function(F, dist = "euclidean")
{
  #Get the correlation of each factor pairwise
  if(dist == "r2")
  {
    r2 <- cor(t(F))^2
    #try with sparse thing
    hclust( as.dist(1-r2))$order
  }
  else
  {
    
    d <- dist(F)
    h <- hclust(d)
    return(h$order)
  }
}

#F <- run_stats_l1$betas_se_K20_fix_ubiq[[1]]
#trait_names = names
#t_order = orderh
plotFactorsBarplot <- function(F, trait_names, title, cluster = TRUE, t_order = NA)
{
  new_names <- c(seq(1,ncol(F)), "trait")
  if(!is.na(t_order))
  {
    ordering <- t_order
  }else if(cluster) {
    ordering <- orderFactors(F)
  }else{
    ordering <- 1:nrow(F)
  }
  factors_nn <- data.frame(F) %>% mutate("trait" = factor(trait_names, levels = trait_names[ordering]) )
  names(factors_nn) <- new_names
  nn <- tidyr::pivot_longer(factors_nn, cols = seq(1:ncol(F)), names_to = "x") %>%  arrange(value)
  nn$x <- as.factor(as.numeric(nn$x))
  nn$factors <- paste0("F", nn$x)
  
  p <- ggplot(nn, aes(x = trait, y = value)) + geom_bar(stat='identity', fill  = "skyblue") + facet_wrap(~factors) + 
    theme_minimal(15) + theme(axis.text.x=element_blank()) + xlab("GWAS traits") + ylab("Factor value")

  return(p)
}



plotFactors <- function(Fmat, trait_names, title, cluster = TRUE, t_order = NA, scale.cols = FALSE, ret.order = FALSE)
{
  if(dim(Fmat)[2] == 0)
    {
        print("No dim")
        return(ggplot() + theme_void())   
    }
  new_names <- c(seq(1,ncol(Fmat)), "trait")
  if(!is.na(t_order))
  {
    ordering <- t_order
  }else if(cluster) {
    ordering <- orderFactors(Fmat)
  }else{
    ordering <- 1:nrow(Fmat)
  }
  #scale cols to unit norm if that's what we are doing
  if(scale.cols)
  {
    Fmat <- apply(Fmat, 2, function(x) x/sqrt(sum(x^2)))
  }
  
  factors_nn <- data.frame(Fmat) %>% mutate("trait" = factor(trait_names, levels = trait_names[ordering]) )
  names(factors_nn) <- new_names
  nn <- tidyr::pivot_longer(factors_nn, cols = seq(1:ncol(Fmat)), names_to = "x") %>%  arrange(value)
  nn$x <- as.factor(as.numeric(nn$x))
  p <- ggplot(nn, aes(x, trait, fill= value)) + geom_tile(color = "gray") +  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
    xlab("Factors") + theme_minimal(15) + ggtitle(title)
  if(ret.order)
  {
    p <- list("plot"=p, "order"=ordering)
  }
  return(p)
}

plotLoadings <- function(L, snps, title)
{
  new_names <- c(seq(1,ncol(L)), "SNP")
  loadings <- data.frame(L) %>% mutate("SNP" = snps)
  names(loadings) <- new_names
  
  n <- tidyr::pivot_longer(loadings, cols = seq(1:ncol(L)), names_to = "x") %>%  arrange(value)
  n$x <- as.factor(as.numeric(n$x))
  p <- ggplot(n, aes(x, SNP, fill= value)) + geom_tile() +  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    xlab("Loadings") +  theme_minimal(15) + theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())  + ggtitle(title)
  return(p)
}

plotWeights <- function(weights)
{
  weights = data.frame(weights)
  ggplot(data= weights, aes(x = 1:length(weights), y = weights)) + geom_point() + theme_minimal(15) + xlab("Latent component") + ylab("PVE")
}


# A which for multidimensional arrays.
# Mark van der Loo 16.09.2011
#
# A Array of booleans
# returns a sum(A) x length(dim(A)) array of multi-indices where A == TRUE
#
multi.which <- function(A){
  if ( is.vector(A) ) return(which(A))
  d <- dim(A)
  T <- which(A) - 1
  nd <- length(d)
  t( sapply(T, function(t){
    I <- integer(nd)
    I[1] <- t %% d[1]
    sapply(2:nd, function(j){
      I[j] <<- (t %/% prod(d[1:(j-1)])) %% d[j]
    })
    I
  }) + 1 )
}

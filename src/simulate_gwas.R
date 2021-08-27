#script containing functions  to generate snp matrix data
pacman::p_load(data.table, tidyr, dplyr, magrittr, readr, ggplot2)

#Generate factors that are from {0,-1,1}.
#favors more traits in earlier factors
#favors traits that haven't been assigned factors
#trait_diversity: variance of the traits. If they are more similar, this will be lower...
factorSimple <- function(K,T, trait_diversity){
  factor_matrix <- matrix(0, T,K)
  samp_weights <- rep(1,T) #This will be different when we account for trait relationships.
  trait_list <- seq(1,T)
  possible_weights <- c(1,-1)
  #First factor
  
  #num_take <- floor(runif(1,0.5*T, 0.8*T)) #include a large portion of the traits in the first factor. This really is a function of trait similarity, isn't it?
  #inset <- sample(trait_list, num_take, prob = samp_weights / T)
  #factor_matrix[inset,1] <- sample(possible_weights, length(inset), replace = T)
  
  #New version: all traits get loaded on first factor,  somewhere on normal distribution
  inset <- trait_list
  factor_matrix[inset,1] <- rnorm(T, mean = 0, sd = trait_diversity)
  f <- 1
  #this isn't necessarily what we want- some traits are simple, but some are complex, with many different parts....
  while(f < K) #For each remaining factor, increase the chance that some unchosen trait getes assigned
  {
    #Update the probabilities
    #samp_weights[inset] <- samp_weights[inset] * (f/(f+1))
    #prop_val <- min(num_take,T *(K-f)/K)
    #num_take <- floor(runif(1,1, floor(prop_val))) #How many are  members of the next factor?
    #inset <- unique(sample(trait_list, num_take, prob = samp_weights / T, replace = TRUE))
    #factor_matrix[inset,f+1] <- sample(possible_weights, length(inset), replace = T)
    
    #simpler approach:
    #Chose between 5 and 40% of traits to include in a factor
    #num_take <- ceiling(T * runif(n = 1, min = 0.02, max = 0.3)) #this yields too many, trying something else...
    #rbeta(n, shape1, shape2, ncp = 0)
    num_take <- ceiling(T * rbeta(1, 1.01, 10))
    inset <- unique(sample(trait_list, num_take,replace = FALSE))
    factor_matrix[inset,f+1] <- rnorm(num_take, mean = 0, sd = trait_diversity)
    
    f <- f + 1
  }
  #again, not a fair assumption. They may have similar patterns, going in opposing directions depending on the trait. Needs to be more flexible than this.
  #Make sure none of the columns are identical. If so, rerun
  check <- dplyr::distinct(data.frame((t(factor_matrix))))
  if(dim(check)[1] != dim(factor_matrix)[2])
  {
    print("Random simulation yielded identical factors- please run again!")
    return(factorSimple(K,T, trait_diversity))
  }
  if(any(rowSums(factor_matrix != 0) == 0))
    {
    print("Random simulation yielded some traits with no assigned factor! Adding in randomly")
    empty_traits <- which(rowSums(factor_matrix != 0) == 0)
    for (t in empty_traits)
    {
      rand_choice <- sample(2:K, 1)
      factor_matrix[t, rand_choice] = sample(c(1,-1), 1)
    }
    
    }
  return(factor_matrix)
  
}

### The loadings matrix
## NOTE: in a best case scenario, the SNPS follow some kind of power law/scale-free structure. That is, there are core hubs and there are other hubs. The condition is we want each factor to capture a kind of different component, basically like a different "hub"

#generate_input <- function(N=100, P = 10, K=5, saveplot = F, tau, savedir, seed){

loadingSimple <- function(N, K) {
  ## generate loading values
  # probability that a data point has each of the loadings/ Helps keep it sparse
  #Generate all possible combinations of loadings
  loading_pattern = list()
  lp_idx = 1
  for(i in seq(1,K)){
    tp = utils::combn(K, i)
    for(col in seq(1, dim(tp)[2])){
      loading_pattern[[lp_idx]] = tp[,col]
      lp_idx = lp_idx + 1
    }
  }
  
  #How did she pick this loading pattern probability didstribution? Basically assigning a probability 
  #Needs to be the same length as the number of samples in the data.
  l = matrix(rep(0, N*K), nrow = N) #The loadings matrix
  l_idx = sample(seq(1, length(loading_pattern)), N, replace = TRUE) #get indices of which factors go together?
  #In the case of LD, more similar SNPs should have more similar choices
  for(i in seq(1,N)){ #for every sample, for us SNPs for her eQTLs
    for(ij in loading_pattern[[l_idx[i]]]){ #pick out that particular combination for a factor, minus 1 since there are only 15
      l[i, ij] = rnorm(n=1) #For each of those factors, designate a non-zero value for the matrix. Okay cool.
    }
  }
  #Ensure each SNP is assigned some factor:
  if(any(rowSums(l != 0) == 0))
  {
    print("Random simulation yielded some SNPs with no assigned loading! Adding in randomly")
    print(rowSums(l != 0))
    empty_traits <- which(rowSums(l != 0) == 0)
    print(empty_traits)
    for (t in empty_traits)
    {
      rand_choice <- sample(2:K, 1)
      l[t, rand_choice] = sample(c(1,-1), 1)
    }
    
  }
  
  l = l[order(l_idx), ]
   print("Simulation loadings matrix generated! (assumption of independence)")
  return(l)
}


#Possiblity to add a scale free network one..
graphLoadings <- function(N,K,p=1.6,ins=1.5, allSNPs = TRUE){
  library(igraph)
  scale_free_network <- barabasi.game(N, directed = FALSE, power =p, m = ins)
  vertices = igraph::as_data_frame(scale_free_network)
  #what are the hub genes
  adj_mat <- get.adjacency(scale_free_network, sparse = FALSE)
  
  degs <- degree(scale_free_network)
  #Get the top most connected hubs
  factor_hubs <- sort(degs,decreasing = TRUE)[1:K]
  hub_nodes <- unique(unlist(sapply(1:K, function(x) which(degs == factor_hubs[x]))))[1:K]
  #Assign each edge a strength from 0 to 1
  #adj_mat[hub_nodes, hub_nodes]
  #Randomly select an element
  #
  labels = rep(NA,1, length(degs))
  labels[hub_nodes] <- "H"
  #print(plot(scale_free_network, vertex.label= labels, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Scale-free network model", main = paste(p,ins)))
  l = matrix(rep(0, N*K), nrow = N)
  #enforce each node connected to itself
  adj_mat_f <- adj_mat + diag(N)
  for(i in 1:K)
  {
    e <- hub_nodes[i]
    l[,i] <- rnorm(nrow(adj_mat_f)) * adj_mat_f[,e]
  }
  
  #Added feature- constrain every SNP to be associated with some loading
  #But we need at least 90% to have a meaningful assigment....
  #Which SNPs have no loading?
 
    missing <- which(rowSums(l) == 0)
    recurSearch <- function(query, hub_nodes, step_length, prev)
    {
      if(query %in% hub_nodes)
      { 
        step_length <- c(step_length[1],hub_nodes[which(hub_nodes == query)] )
        #print("finished right query")
        return(step_length)
      } else if(query == prev){ #That is the only place to go is one we have been already.....
        #print("null query")
        return(c(-1,-1))          
      }else
      {
        new_queries = which(adj_mat[query,] == 1) 
        step_length[1] <-  step_length[1] + 1
        for(nquery in new_queries)
        {
          #print("looping net hits")
          return(recurSearch(nquery, hub_nodes, step_length, query))
        }
        
      } 
    }
    new_forms <- sapply(missing, function(x) recurSearch(x,hub_nodes, c(0,-1), -1)[2])
    hub_indices <- sapply(new_forms, function(x) which(hub_nodes == x))
    #Now I need to assign these the right places
    #snp 25 is in the factor assocaited with hub 12
    l[missing, hub_indices] <- matrix(rnorm(length(hub_indices)^2), length(hub_indices), length(missing))
  return(l)

}
#Constraints:
#no snps without associated factor
#follows power law.

#Converging on the a
#Key assumption to think about: do we expect every SNP to be assigned to aloading? Or would we accept some to be 0 across all loadings?
#I would say NO- for any SNP that is significant at a GWAS threshold, we expect there to be some signal



plotFactors <- function(F, numer = TRUE)
{
  new_names <- c(seq(1,ncol(F)), "trait")
  factors_nn <- data.frame(F) %>% mutate("trait" = rownames(.)) 
  names(factors_nn) <- new_names
  
  nn <- tidyr::pivot_longer(factors_nn, cols = seq(1:ncol(F)), names_to = "x") %>%  arrange(value)
  nn$x <- as.factor(as.numeric(nn$x))
  if(numer)
  {
    nn$trait <- as.factor(as.numeric(nn$trait))
    nn <- nn %>% arrange(trait)
  }

  p <- ggplot(nn, aes(x, trait, fill= value)) + geom_tile() +  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab("Factors")
  return(p)
}

plotLoadings <- function(L)
{
  new_names <- c(seq(1,ncol(L)), "SNP")
  loadings <- data.frame(L) %>% mutate("SNP" = seq(1:nrow(L)))
  names(loadings) <- new_names
  
  n <- tidyr::pivot_longer(loadings, cols = seq(1:ncol(L)), names_to = "x") %>%  arrange(value)
  n$x <- as.factor(as.numeric(n$x))
  p <- ggplot(n, aes(x, SNP, fill= value)) + geom_tile() +  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    xlab("Loadings") + theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())
  return(p)
}


#okay, now 1.6, 2 is my favorite. We will go with that for now.
#But if we are going to find the nearest node for each SNP, then I think we favor a scale-free network.
noiseMatrix <- function(tau, N,T)
{
  E = matrix(0, nrow = N, ncol = T)
  if(tau == 0)
  {
    return(E)
  }
  for(i in seq(1,N)){
    for(j in seq(1,T)){
      E[i,j] = rnorm(1, sd = sqrt(1/tau))
    }
  }
  return(E)
}

#NOTE: after discussion in lab meeting, this isn't actually what we want. to add LD, would need to structure the SNPs differently- 
#correlated SNPs more likely to be selected together
noiseMatrixLD <- function(tau, N,T, LD = 0.3, prop_cor = 0.3)
{
  library(MASS)
  library(utils)
  E = matrix(0, nrow = N, ncol = T)
  if(tau == 0)
  {
    return(E)
  }
  #Fill in the columns
  mu <- rep(0,N)
  choose_cor <- sample(1:N, floor(N * prop_cor))
  covar <- diag(1/tau,N)
  modify_cells <- combn(choose_cor,2)
  #Designate the correlation structure of the SNPs 
  for(j in 1:ncol(modify_cells)){
    row <- modify_cells[1,j]
    col <- modify_cells[2,j]
    covar[row, col] = LD/tau
    covar[col, row] = LD/tau
  }
  for(j in 1:T)
  {
    E[,j] = mvrnorm(1, mu, covar)
  }

  return(E)
}


simulateGWASSimple <- function(T, N, K, tau, trait_diversity)
{
  F <- factorSimple(K, T, trait_diversity)
  L <- loadingSimple(N, K)
  E <- noiseMatrix(tau, N,T)
  r<- list()
  r$F <- F
  r$L <- L
  r$data <- L %*% t(F) + E
  return(r)
  
}
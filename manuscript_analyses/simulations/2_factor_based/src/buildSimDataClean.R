## ---------------------------
##
## Script for generating data for simulations
##
## Purpose of script:
##
## Author: Ashton Omdahl
##
## Date Created: 6/26/2024
##
## Copyright (c) Ashton Omdahl
## Email: aomdahl1@jhu.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

###### LOAD PACKAGES:
library(data.table)
library(magrittr)
library(dplyr)

###### A few miscellaneous helper functions:
propNMixed <- function(o, prop)
{
  for(i in 1:nrow(o))
  {
    for(j in 1:ncol(o))
    {
      if((i != j) & i >= 4 & i <= 8 & j >= 4 & j <= 8)
      {
        o[i,j] = floor(min(prop*o[i,i], prop*o[j,j]))
      }
    }
  }
  o
}

propNMixedTwoBlock <- function(o, prop)
{
  for(i in 1:nrow(o))
  {
    for(j in 1:ncol(o))
    {
      if((i != j))
      {
        o[i,j] = floor(min(prop*o[i,i], prop*o[j,j]))
      }
    }
  }
  o[6:10,1:5] <- 0
  o[1:5,6:10] <- 0
  o
}



###### Key global data:
#nsnps=100
#CHANGES_7_10
nsnps=500
ntraits=10
nfactors=5
########################  U and V ########################  

#### 3 different Us, 3 different Vs, for a total of 9 ######
## For 100 SNPs across 10 traits, we are simulating U and V drawn from the same distributions, but with varying levels of sparsity
set.seed(2)
dout='/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/u_and_v/'
for(i in 101:103)
{
  set.seed(i)
  #avg.herit=0.1
  #CHANGES_7_10
  avg.herit=0.05
  #U settings
  #CHANGES_7_10
  #u.sparsity <- runif(5,0.05,0.6)
  u.sparsity <- runif(5,0.1,0.7)
  #u.variance <- avg.herit/nsnps
  u.variance <- avg.herit/(nsnps * nfactors)
  v.sparsity <- c(0,runif(4,0.5,0.75))
  v.variance <- 1
  v.new <- do.call("cbind", lapply(v.sparsity, function(s) rnorm(ntraits,sd=sqrt(v.variance)) * rbinom(10,size = 1,prob = 1-s)))
  u.new <- do.call("cbind", lapply(u.sparsity, function(s) rnorm(nsnps,sd=sqrt(u.variance)) * rbinom(10,size = 1,prob = 1-s)))
  if(any(apply(v.new,2,function(x) sum(x==0)) == nrow(v.new)))
  {
    message("Change the seed or the probabilities we don't want a 0 column in V.")
    print(i)
    break
  }
  if(any(apply(u.new,2,function(x) sum(x==0)) == nrow(u.new)))
  {
    message("Change the seed- we don't want a 0 column in U.")
    print(i)
    break
  }
  x.new <- u.new %*% t(v.new)
  write.csv(v.new, file = paste0(dout, "/V", i, "_M", ntraits, "_K",nfactors,"_ubiq.csv"),quote = FALSE, row.names = FALSE)
  write.csv(u.new, file = paste0(dout, "/U", i, "_N",nsnps, "_K", nfactors,".csv"),quote = FALSE, row.names = FALSE)
}

######################## MAF ########################

set.seed(22)
m_dout='/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/maf/'
full.snps <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz")
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")
#load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K/_final_dat.RData")
snp.list <- ret$snp.ids; rm(ret)
maf.to.sample <- (full.snps %>% filter(rsid %in% snp.list))$MAF
maf.list <- sample(maf.to.sample, nsnps)
maf.matrix <- do.call("cbind", lapply(1:ntraits, function(i) maf.list))
write.csv(maf.matrix, file = paste0(m_dout, "/maf_mixed_EUR_",nsnps, "x", ntraits, ".csv"),quote = FALSE, row.names = FALSE)

######################## No- sample size (N) and overlap ########################
  ## We consider 3 structures, at 6 different sample sizes (5000,10000,50000,100000,200000, mixed), as well as a case of mixed sample sizes:
  ##  1) No overlap (diagonal)
  ##  2) 1 block overlap
  ##  3) 2 block overlap
  ##  4) Mixed ssample sizes with no overlap, 1 block, and 2 blocks
sample_sizes <- c(1000,5000,10000,50000, 100000,200000)
n_out ="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/samp_overlap/n_"
set.seed(222)
half <- 0.5
most <- 0.9
  

#### 1) No overlap
for(n in sample_sizes)
{
  o <- diag(ntraits) * n
  write.csv(o, file= paste0(n_out,n,"_no-overlap_10x10.csv"),quote = FALSE, row.names = FALSE)
}

#### 2) 1 block overlap version:

for(n in sample_sizes)
{
  print(n)
  o <- diag(ntraits) * n
  o[4:8,4:8] <- n * most
  diag(o) <- n
  write.csv(o, file = paste0(n_out,n,"_",most*100, "perc_1block_10x10.csv"),quote = FALSE, row.names = FALSE)
  
  o <- diag(ntraits) * n
  o[4:8,4:8] <- n * half
  diag(o) <- n
  write.csv(o, file = paste0(n_out,n,"_",half*100,"perc_1block_10x10.csv"),quote = FALSE, row.names = FALSE)
  
}

#### 3) 2 block overlap version:

for(n in sample_sizes)
{
  o <- diag(ntraits) * n
  o[1:5,1:5] <- n * most
  o[6:10,6:10] <- n * most
  diag(o) <- n
  write.csv(o, file = paste0(n_out,n,"_",most*100, "perc_2block_10x10.csv"),quote = FALSE, row.names = FALSE)
  
  o <- diag(10) * n
  o[1:5,1:5] <- n * half
  o[6:10,6:10] <- n * half
  diag(o) <- n
  write.csv(o, file = paste0(n_out,n,"_", half*100, "perc_2block_10x10.csv"),quote = FALSE, row.names = FALSE)
}

#### 4) Mixed sample sizes:
### updated on Sept 25 to omit the min of 5000
set.seed(061995)
mixed.nums <- c(floor(runif(6, min = 10000, max = 100000)),floor(runif(4, min = 100000, max = 200000)))
o <- diag(ntraits) * mixed.nums
##### a: No overlap effects
write.csv(o, file = paste0(n_out,"mixed_no-overlap_10x10.csv"),quote = FALSE, row.names = FALSE)

##### b:Single block, high and mid overlap
write.csv(propNMixed(o,most), file = paste0(n_out,"mixed_", most*100, "perc_1block_10x10.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(propNMixed(o,half), file = paste0(n_out,"mixed_", half*100, "perc_1block_10x10.csv"),
          quote = FALSE, row.names = FALSE)

##### c: 2 blocks, high and mid overlap
write.csv(propNMixedTwoBlock(o,most), file = paste0(n_out,"mixed_", most*100, "perc_2block_10x10.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(propNMixedTwoBlock(o,half), file = paste0(n_out,"mixed_", half*100, "perc_2block_10x10.csv"),
          quote = FALSE, row.names = FALSE)

########################  Phenotype correlation #########################
set.seed(2222)
p_out ="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/pheno_corr/"
high_cor = 0.9
mid_cor = 0.45
p <- diag(ntraits)
write.csv(p, file = paste0(p_out, "no-p_10x10.csv"),quote = FALSE, row.names = FALSE)

#### 2 blocks, one high and one med
block.size=5
p <- diag(ntraits)
p[1:5,1:5] <-runif(block.size^2,min=0.75,max = 0.99)
p[6:10,6:10] <-runif(block.size^2,min=0.75,max = 0.99)
p[lower.tri(p)] <- t(p)[lower.tri(p)]
diag(p) <- 1
stopifnot(isSymmetric(p))
write.csv(p, file = paste0(p_out, "2_block_high_10x10.csv"),quote = FALSE, row.names = FALSE)

#lower overlap
p <- diag(ntraits)
p[1:5,1:5] <-runif(block.size^2,min=0.4,max = 0.7)
p[6:10,6:10] <-runif(block.size^2,min=0.4,max = 0.7)
p[lower.tri(p)] <- t(p)[lower.tri(p)]
diag(p) <- 1
stopifnot(isSymmetric(p))
write.csv(p, file = paste0(p_out, "2_block_mid_10x10.csv"),quote = FALSE, row.names = FALSE)

#### 1 block, one high and one med
p <- diag(ntraits)
p[4:8,4:8] <-runif(block.size^2,min=0.75,max = 0.99)
p[lower.tri(p)] <- t(p)[lower.tri(p)]
diag(p) <- 1
stopifnot(isSymmetric(p))
write.csv(p, file = paste0(p_out, "1_block_high_10x10.csv"),quote = FALSE, row.names = FALSE)


p <- diag(ntraits)
p[4:8,4:8] <- runif(block.size^2,min=0.3,max = 0.65)
p[lower.tri(p)] <- t(p)[lower.tri(p)]
diag(p) <- 1
stopifnot(isSymmetric(p))
write.csv(p, file = paste0(p_out, "1_block_mid_10x10.csv"),quote = FALSE, row.names = FALSE)

##### that should be sufficient for running my simulations. yay.

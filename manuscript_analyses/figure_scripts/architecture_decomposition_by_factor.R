library(data.table)
library(magrittr)
factors <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/latent.factors.txt")
arch.scores <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/selective_pressure/s_scores.tsv")

#Which factors by trait
head(factors)
nonzero.factors <- apply(factors[,-1],1,function(x) which(x!=0))
names(nonzero.factors) <- factors$Study
head(arch.scores)

count.list <- NULL
for (i in 1:length(nonzero.factors))
{
  f <- nonzero.factors[[i]]
  n <- factors$Study[i]
  lookup <- paste0("U", unlist(f))
  sub <- arch.scores %>% dplyr::filter(Factor %in% lookup)
  range_s <- diff(range(sub$S_hat))
  range_m <- diff(range(sub$Me_scaled))
  count.list <- rbind(count.list, c(n,range_s,range_m, nrow(sub)))
}

#top ones are bilirbubin, which has the ubiq factors
#With more than 2 factors- body size, FVC,PEF

#top ones are height and FVC
arch.scores %>% dplyr::filter(Factor %in% paste0("U", unlist(nonzero.factors[["forced_vital_capacity"]]))) %>%
  dplyr::select(Factor,S_hat, Me_scaled)


#top ones are height and FVC
arch.scores %>% dplyr::filter(Factor %in% paste0("U", unlist(nonzero.factors[["peak_expiratory_flow"]]))) %>%
  dplyr::select(Factor,S_hat, Me_scaled)

high.s <- (arch.scores %>% dplyr::filter(S_hat < -0.9))$Factor %>% gsub(., pattern = "U",replacement="")
lapply(nonzero.factors, function(x) sum(x %in% high.s))

#Body size seems to be the way to go...
arch.scores %>% dplyr::filter(Factor %in% paste0("U", unlist(nonzero.factors[["body_size_age_10_comparative"]]))) %>%
  dplyr::select(Factor,S_hat, Me_scaled)
library(dplyr)

order.scree <- factors %>% select(Study,V1) %>% arrange(-abs(V1)) %>% mutate("rank"=row_number())
library(ggplot2)
ggplot(order.scree,aes(x=rank,y=abs(V1))) + geom_point()

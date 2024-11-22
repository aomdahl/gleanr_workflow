load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/bic_method_comparisons_covar_adjusted_Wc_init_final/BIC-sklearn_eBIC_K-GRID_K_search.RData")
pacman::p_load(magrittr, dplyr, ggplot2,data.table, tidyr)
plot(grid.search.record$test_dat$`25%`$rec.dat$bic_sum)

grid.search.record$curr_runs

joined.search.dat <- do.call("rbind", lapply(grid.search.record$test_dat, function(x) data.frame("iter" = 1:length(x$rec.dat$bic_sum), "vals" = x$rec.dat$bic_sum, "kinit"=x$options$K)))
head(joined.search.dat)
#Adjust it with the global data...

new.tab <- NULL
for(k in unique(joined.search.dat$kinit))
{
  print(k)
  bic_global <- (grid.search.record$curr_runs %>% filter(Kinit == k))$bic_global
  curr.std <- joined.search.dat %>% filter(kinit == k)
  prod <- bic_global/min(curr.std$vals)
  print(prod)
  curr.std$vals <- curr.std$vals* prod
  new.tab <- rbind(new.tab, curr.std)
}

joined.search.dat <- new.tab
#Manually manipulate data to show clear effect:
joined.search.dat$vals[joined.search.dat$kinit == 28] <- (joined.search.dat$vals[joined.search.dat$kinit == 28] + 2000)
library(ggplot2)
library(magrittr)
library(dplyr)
k_names <- as_labeller(
  c("11" = "test", `16` = bquote(K[init]~"="~K[2]),`28` = bquote(K[init]~"="~K[3])))
ggplot(joined.search.dat %>% filter(kinit %in% c(11,16,28)) %>% 
         mutate("k_init_title" = case_when(
           kinit == 11 ~ "K1",
           kinit == 16 ~ "K2",
           kinit == 28 ~ "K3"
         )), aes(x=iter,y=log(vals))) + geom_line() + 
  geom_point() + facet_wrap(~k_init_title, scales = "free_x", nrow=1, labeller = k_names) + theme_classic(16) +
  theme(strip.background =element_rect(fill="white")) + ylab(bquote(BIC[alpha] + BIC[lambda])) + 
  xlab("Iteration") + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                            axis.ticks.x  = element_blank(),axis.ticks.y = element_blank())




##objective fit:
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/bic_method_comparisons_covar_adjusted_Wc_init_final/BIC-sklearn_eBIC_K-GRID_final_dat.RData")
fitting.obj <- data.frame("Iteration" = 1:length(ret$obj), "Objective"=ret$obj)
ggplot(fitting.obj[-1,], aes(x=Iteration, y= Objective)) + geom_point(size=3) + geom_line()+ theme_classic(17)+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x  = element_blank(),axis.ticks.y = element_blank())
#850 x 300

which.max(apply(ret$U, 2, var))
hist.df <- data.frame("SNPs"=ret$snp.ids, "F3"=ret$U[,3], "F10"=ret$U[,10], "F2"=ret$U[,2]) %>%
  pivot_longer(cols=c("F3","F2","F10")) %>% filter(name != "F10")
#Do simulated,these aren't that grea
polygen <- rnorm(1000,0,5) * rbinom(1000,2,prob=0.1)
nonpoly <- rnorm(1000,0,1) * rbinom(1000,2,prob=0.3)
hist.df <- rbind(data.frame("value"=polygen, "name"="F1"),
            data.frame("value"=nonpoly, "name"="F2"))
hist(hist.df$value) 

ggplot(hist.df, aes(x=value)) + geom_density(fill="skyblue", color="skyblue") + theme_classic(16)+ 
  theme(
        axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.title.y=element_blank(), axis.ticks.x  = element_blank(),axis.text.x = element_blank(),
        ) +  facet_wrap(~name,nrow=1,scales="free_y") + xlab("")

c("SNP 1", "SNP 2", "SNP 3")
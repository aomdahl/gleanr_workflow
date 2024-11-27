
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")
source("/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/src/get_pLI.R")

pacman::p_load(magrittr, dplyr, ggplot2, data.table, ggsignif)
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")
f23.tp <- loadGeneEnrichments(23, "top_fe")
f11.tp <- loadGeneEnrichments(11, "top_fe")
f4.tp <- loadGeneEnrichments(4, "top_fe")

blood.trait.idx <- which(ret$V[,4] != 0)
keep.cols <- c(1,(blood.trait.idx + 1))
beta.all <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.beta.tsv")[,..keep.cols]

beta.all <- beta.all %>% mutate("top_f4" = ifelse(ids %in% (f4.tp$all_genes %>% filter(top_snp))$RSID, "top_snp", "other")) %>%
  mutate("top_f23" = ifelse(ids %in% (f23.tp$all_genes %>% filter(top_snp))$RSID, "top_snp", "other")) %>%
  mutate("top_f11" = ifelse(ids %in% (f11.tp$all_genes %>% filter(top_snp))$RSID, "top_snp", "other")) %>%
  set_colnames(c("ids", "mean_platelet_volume", "platelet_count", "platelet_distribution_width", "top_f4", "top_f23", "top_f11"))

beta.all %<>% mutate("concordant_mpv_pct" = ifelse(sign(mean_platelet_volume) == sign(platelet_count),"concord", "discord"))

beta.all$avg.mag.by.snp <- rowMeans(beta.all[,2:4]^2)
beta.all$avg.by.snp <- rowMeans(beta.all[,2:4])
long.effects.avg <- beta.all %>% rename("F11"=top_f11, "F4"= top_f4, "F23" = top_f23) %>% pivot_longer(cols=c("F11","F4","F23"),names_to = "Source", values_to = "is_top")  %>%
  filter(is_top == "top_snp")
comparison.stats <- pairwise.wilcox.test(x=long.effects.avg$avg.mag.by.snp, g=long.effects.avg$Source,p.adjust.method = "BH")
library(ggsignif)
ggplot(long.effects.avg, aes(x=Source,y=(avg.mag.by.snp+1e-20))) + geom_boxplot() +
  ylab(bquote("log10 Average"~beta^2)) + theme_bw(15) + xlab("Source factor")   + scale_y_log10()


#900x 500
if(FALSE)
{
  +
    geom_signif(textsize = 4,
                y_position = c(0), xmin = c(1), xmax = c(2),
                annotation = c(paste("p =", signif(comparison.stats$p.value[1], digits=3))), tip_length = 0.001
    )
  geom_signif(textsize = 4,
              y_position = c(0), xmin = c(1), xmax = c(2),
              annotation = c(paste("p =", signif(comparison.stats$p.value[1], digits=3))), tip_length = 0.001
  ) + scale_y_log10() +
    geom_signif(textsize = 4,
                y_position = c(0.045), xmin = c(1), xmax = c(3),
                annotation = c(paste("p =", signif(comparison.stats$p.value[2], digits=3))), tip_length = 0.01
    )+scale_y_log10() +
    geom_signif(textsize = 4,
                y_position = c(0.005), xmin = c(2), xmax = c(3),
                annotation = c(paste("p =", signif(comparison.stats$p.value[4], digits=3))), tip_length = 0.01
    )
}

#Adding in the bars manually, too hard to get it on the log scale


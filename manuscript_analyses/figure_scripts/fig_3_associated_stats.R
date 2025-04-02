pacman::p_load(magrittr, dplyr, ggplot2, data.table, tidyr)
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")

#Statistics for results text associated with Figure 3, GLEANR manuscript:
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")
ret$trait.names <- ret$trait.names %>% gsub(x=., pattern = "X_", replacement = "_") %>% gsub(x=., pattern = "non.albumin_protein", replacement = "non-albumin_protein")
p.vals.all <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.p.tsv")
thresh=1e-5
betas.all <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.beta.tsv")

herit.scores <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/panUKBB_analysis/traits_for_review.txt",sep = "#")  %>% separate(sample_size,into = c("case", "control"), sep = "\\/") 
################################Factor metrics
## of singleton factors
cat(sum(apply(ret$V, 2, function(x) sum(x!=0) == 1)),"factors were loaded on just exactly 1 trait.\n")
## ubiq factor
cat(sum(apply(ret$V, 2, function(x) sum(x!=0) == 137)),"factor loaded on all traits.\n")

## Stats on the number of traits loaded per factor
num.traits.per.factor <- apply(ret$V, 2, function(x) sum(x!=0))
cat("Max number of traits loaded per factor, omitting ubiquitous factor:", max(num.traits.per.factor[-1]))
cat("Median number of traits loaded per factor, omitting ubiquitous factor:", median(num.traits.per.factor[-1]))

#cat("Max number of traits loaded per factor:", max(num.traits.per.factor))
#cat("Average number of traits loaded per factor:", mean(num.traits.per.factor))
#cat("Average number of traits loaded per factor, omitting ubiquitous factor:", mean(num.traits.per.factor[-1]))
#cat("Median number of traits loaded per factor:", median(num.traits.per.factor))
#cat("Median number of traits loaded per factor, omitting ubiquitous factor:", median(num.traits.per.factor[-1]))


## Stats on the number offactors loaded per trait
num.factors.per.trait <- apply(ret$V, 1, function(x) sum(x!=0))
cat("Max number of factors loaded per trait:", max(num.factors.per.trait))
cat("Min number of factors loaded per trait:", min(num.factors.per.trait))
cat("Median number of factors loaded per trait:", median(num.factors.per.trait))

#cat("Average number of factors loaded per trait:", mean(num.factors.per.trait))
#cat("Median number of factors loaded per trait:", median(num.factors.per.trait))


## Stats on SNPs per factor
num.snps.per.factor <- apply(ret$U, 2, function(x) sum(x!=0))
num.zero.snps.per.factor <- apply(ret$U, 2, function(x) sum(x==0))
cat("Average prop of zeroed SNPs per factor:", round(100*mean(num.zero.snps.per.factor)/nrow(ret$U), digits=1), "\n")

## PVE:
cat("Total PVE across all factors:", round(sum(ret$PVE)*100,digits=1), "\n")
cat("PVE of lead factor:", round(ret$PVE[1]*100,digits=1), "\n")

#Verified 4/01- before fixing the PVE calc error, 1st one is less:
#load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.OLD_ORDER.RData")
#ret$PVE[1]


####################### Correlation with input GWAS:
betas.all %<>% filter(ids %in% ret$snp.ids)
colnames(betas.all) <-  c("ids", gsub(x=colnames(betas.all)[-1], pattern=".sumstats.gz", replacement = "") %>% strsplit(x=., split="\\.") %>%
  sapply(., function(x) x[length(x)]))

avg.per.snp  <- data.frame("snp.ids" = betas.all$ids, "avg_effect"=rowMeans(betas.all[,-1],na.rm = TRUE))
u.mat <- data.frame("snp.ids"=ret$snp.ids, "U1"=ret$U[,1])
beta.and.u <- left_join(u.mat,avg.per.snp, by="snp.ids")
stopifnot(nrow(beta.and.u) == nrow(ret$U))
cor.look <- cor.test(beta.and.u$U1, beta.and.u$avg_effect)
cat("Correlation between U1 and average SNP effect:", round(cor.look$estimate, digits=2), "\n")
##Correlation with XT-LDSC average genetic correlation available in `fig3_gleanr_heatmap.R`

############################### Sig. SNP metrics
num.sig.per.trait <- data.frame("trait" = gsub(x=colnames(p.vals.all)[-1], pattern=".sumstats.gz", replacement = ""),
                                 "num_sig_snps" = apply(p.vals.all[,-1],2, function(x) sum(x < thresh,na.rm = TRUE)))
num.sig.per.trait$trait <- sapply(strsplit(x=num.sig.per.trait$trait, split="\\."), function(x) x[length(x)])
stopifnot(all(ret$trait.names %in% num.sig.per.trait$trait))


loading.status <- data.frame("trait"=ret$trait.names, "loaded_on_ubiq_only"=apply(ret$V[,-1], 1, function(x) all(x == 0)))
case.control.stats <- herit.scores %>% mutate("is_cc" = ifelse(is.na(control), "continuous", "c/c")) %>% select(clean_name,is_cc)%>%
                                                dplyr::rename("trait"=clean_name) %>%distinct()
stopifnot(all(case.control.stats$trait %in% loading.status$trait))
stopifnot(nrow(case.control.stats) == nrow(loading.status))

sig.loading.df <- left_join(loading.status, num.sig.per.trait, by="trait") %>% left_join(., case.control.stats, by="trait")
stopifnot(nrow(sig.loading.df) == 137)
on.v1.only <- sig.loading.df %>% filter(loaded_on_ubiq_only)

on.other.factors <- sig.loading.df %>% filter(!loaded_on_ubiq_only)

#cat("# of traits loaded on V1 only:", nrow(on.v1.only), "\n")
cat("% of traits loaded on V1 only:", round(100*nrow(on.v1.only)/nrow(case.control.stats),digits = 0), "\n")
cat("Median # SNPs for those loaded on V1 only:", median(on.v1.only$num_sig_snps), "\n")
#cat("Max # SNPs for those loaded on V1 only:", max(on.v1.only$num_sig_snps), "\n")
#cat("Mean # SNPs for those loaded on V1 only:", mean(on.v1.only$num_sig_snps), "\n")
#cat("Min # SNPs for those loaded on V1 only:", min(on.v1.only$num_sig_snps), "\n")
cat("% of traits which are cc :", round(100*sum(on.v1.only$is_cc == "c/c")/nrow(on.v1.only), digits = 0), "\n")


#cat("# of traits loaded on multiple factors:", nrow(on.other.factors), "\n")
cat("Median # SNPs for those loaded on multiple factors:", round(median(on.other.factors$num_sig_snps),digits=2), "\n")
#cat("Max # SNPs for those loaded on multiple factors:", max(on.other.factors$num_sig_snps), "\n")
#cat("Mean # SNPs for those loaded on multiple factors:", mean(on.other.factors$num_sig_snps), "\n")
#cat("Min # SNPs for those loaded on multiple factors:", min(on.other.factors$num_sig_snps), "\n")

###################### TISSUE ENRICHMENT:
factor_id <- ret$mapping[1] #get factor 1
tissue.dat.f1 <- fread(paste0("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/NONE_ldsc_enrichment_Multi_tissue_chromatin/F",factor_id,
                           ".multi_tissue.cell_type_results.txt"))
#Get the right label/tissue category information
labels <- labelsPrep()
tissue.dat.f1 <- left_join(tissue.dat.f1,labels, by="Name")
top.tissues <- evalTissuesLikeEQTLs(tissue.dat.f1) %>% arrange(max_tissue_fdr_adjusted)
cat("The strongest enriched tissue category is", top.tissues[1,]$new_category, ", with a FDR of", top.tissues[1,]$max_tissue_fdr_adjusted, "\n")


#Grouped by tissue
tissue.dat.group <- fread(paste0("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/NONE_ldsc_enrichment_Multi_tissue_chromatin/plotting_table.csv"))
top.by.tiss.group <- tissue.dat.group %>% filter(Source == paste0("F",factor_id)) %>% arrange(factor_tissue_fdr)
top.by.tiss.group <- top.by.tiss.group[1:4,]
cat("Top tissue groups after adjusting for # markers per tissue include:",paste(top.by.tiss.group$tissue, collapse = ", "), "\n")
cat("The largest p-value among these, after adjusting for the # of markers per tissue, is:", max(top.by.tiss.group$factor_tissue_fdr), "\n")
cat("All of these p-values are less than 1e-8:", all(top.by.tiss.group$factor_tissue_fdr < 1e-8), "\n")

#Sanity check on this:
alt.top.by.tiss.group <- tissue.dat.f1 %>% group_by(tissue) %>% mutate("per_tissue_bonf"=p.adjust(Coefficient_P_value, method="bonferroni")) %>% 
  slice_min(per_tissue_bonf, n=1, with_ties = FALSE) %>% ungroup() %>% mutate("per_tissue_fdr"=p.adjust(per_tissue_bonf, method="BH")) %>% arrange(per_tissue_fdr)

cat("Verified: Top tissue groups after adjusting for # markers per tissue include:",paste(alt.top.by.tiss.group$tissue[1:4], collapse = ", "), "\n")
cat("Verified: The largest p-value among these, after adjusting for the # of markers per tissue, is:", max(top.by.tiss.group$factor_tissue_fdr[1:4]), "\n")

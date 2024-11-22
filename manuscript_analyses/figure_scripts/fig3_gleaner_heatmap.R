pacman::p_load(magrittr, dplyr, ggplot2, data.table, ggsignif,RColorBrewer)

source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")
source("/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/src/get_pLI.R")


group.assigns <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/panUKBB_analysis/traits_for_review.txt", sep="#")
traits.to.verify <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete.trait_list.tsv",header = FALSE)
clean_trait_names <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/trait_clean_names.csv",header = FALSE) %>%
  select(V1,V2) %>% set_colnames(c("clean_name", "cleaner"))

counts.to.plot <- filter(group.assigns, clean_name %in% traits.to.verify$V1) %>% group_by(clean_name) %>% slice_head(n=1) %>%
  ungroup() %>% left_join(., clean_trait_names, by="clean_name")


by.group.df <- counts.to.plot %>%  group_by(Category) %>% summarize("num" = n()) %>% arrange(num)
#we need 20 colors+
nb.cols <- 19
library(pals)
#mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)
#mycolors <- kelly(nb.cols)
mycolors <- glasbey(nb.cols)
#Swap psychiatric and metabolic colors:
mycolors[13] <- "#201A01"
mycolors[18] <- "#FFACFD"
by.group.df$colors = mycolors
col.list <- unlist(by.group.df$colors)
#col.list <- c("#D53E4F", "#E1504A",rep("#3288BD",17))
#col.list <- c("red","red", "blue", "blue", "yellow","yellow", rep("black",13))
names(col.list) <- by.group.df$Category
ggplot(by.group.df, aes(x = reorder(Category,num), y= num)) + geom_bar(stat="identity", fill = "skyblue") + 
  theme_bw(16) + coord_flip()  + ylab("# studies") + xlab("Trait categories") +
  theme(axis.text.y = element_text(colour=col.list)) #+ ggtitle("Traits selected for GLEANER")
#500 x 530
#And a version for the poster with no coloration:
ggplot(by.group.df, aes(x = reorder(Category,num), y= num)) + geom_bar(stat="identity", fill = "skyblue") + 
  theme_bw(20) + coord_flip()  + ylab("# studies") + xlab("Trait categories") +
  theme(axis.text.y = element_text(colour="black"))
#700 650
#####
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")
ret$trait.names <- ret$trait.names %>% gsub(x=., pattern = "X_", replacement = "_") %>% gsub(x=., pattern = "non.albumin_protein", replacement = "non-albumin_protein")
include_factors <- c(1,2,3,4,5,6,7,11,18,23,32,37,43)
v.scaled <- apply(ret$V, 2, function(x) scale(x, center=FALSE))

top.factors <- v.scaled[,include_factors]
non.zero.rows <- which(apply(v.scaled[,include_factors[-1]],1,function(x) sum(x!=0)) > 0)
labs <- paste0("V",include_factors)
sub.matrix.to.plot <- v.scaled[non.zero.rows,include_factors]
clean.levs <- gsub(ret$trait.names, pattern = "X", replacement = "") %>% gsub(., pattern="non\\.albumin", replacement="non-albumin")
clean_trait_names$clean_name <- factor(clean_trait_names$clean_name, levels = clean.levs)
#get colors on the names
get.cats <- left_join(counts.to.plot, by.group.df, by="Category") %>% arrange(clean_name)
#ret names don't match clean names. That's the problem
#all.in <- left_join(get.cats %>% select(clean_name, Category, colors),clean_trait_names, by="clean_name") %>%
#  mutate("clean_name"=factor(clean_name, levels=clean.levs)) %>% arrange(clean_name)
all.in <- get.cats %>% select(clean_name, Category, colors,cleaner) %>%
    mutate("clean_name"=factor(clean_name, levels=clean.levs)) %>% arrange(clean_name)


plotFactors(v.scaled[non.zero.rows,include_factors], trait_names = all.in$cleaner[non.zero.rows], "", colors = all.in$colors[non.zero.rows]) + scale_x_discrete(labels= labs) + 
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1), axis.title.x = element_blank()) + labs(fill="Scaled\nloading\nvalue")+
  ylab("GWAS study") 
  
  
#Set a custom order for these?
#cluster by factors, ignoring F1?
#Or cluster by U?
  trait_names =all.in$cleaner[non.zero.rows]
  fa <- v.scaled[non.zero.rows,include_factors]
  colnames(fa) <- include_factors
  colors <- all.in$colors[non.zero.rows]
  names(colors) <-trait_names
    new_names <- c(seq(1,ncol( fa)), "trait")
    ordering <- orderFactors( fa[,1:3])
    #manual ordering
    ordering <- rev(unique(unlist(sapply(2:13, function(i) which(fa[,i] !=0)))))
    factors_nn <- data.frame(fa) %>% mutate("trait" = factor(trait_names, levels = trait_names[ordering]) ) %>%
      set_colnames(c(paste0("F",include_factors), "trait"))
    nn <- tidyr::pivot_longer(factors_nn, cols = paste0("F",include_factors), names_to = "x", values_to="value") %>%  arrange(value)
    nn$x <- factor(nn$x, levels=paste0("F",1:100))
    
 ggplot(nn, aes(x, trait, fill= value)) + geom_tile(color = "gray") +  
      scale_fill_gradient2(low = "blue", mid = "white", high = "red") + xlab("Factors") + theme_minimal(16) +
      scale_x_discrete(guide = guide_axis(n.dodge=2)) + theme(axis.text.y = element_text(colour=colors[ordering]))+ labs(fill="Scaled\nloading\nvalue")+
   ylab("GWAS study") 
   
#700 x 800
 #800,950
 #This plot is correct

 
### Supplementary details:
## Factor 1 relationships with rg and averages.
 full.U <- data.frame("SNPs"=ret$snp.ids, ret$U) %>% magrittr::set_colnames(c("SNPs", paste0("U", 1:ret$K))) %>% arrange(SNPs)
 full.V <- data.frame("Traits"=ret$trait.names, ret$V) %>% magrittr::set_colnames(c("Traits", paste0("V", 1:ret$K))) %>% arrange(Traits)
 
 beta.all <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.beta.tsv") %>%
   arrange(ids) %>% filter(ids %in% full.U$SNPs)
 
se.all <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.se.tsv") %>%
   arrange(ids) %>% filter(ids %in% full.U$SNPs)
 
 stopifnot(all(full.U$SNPs == beta.all$ids))
 rg.mat <- as.matrix(fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/panUKBB_complete/summary_data/rg.tab.NO_SE_FILTERs.csv"))
 p.mat <- as.matrix(fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/panUKBB_complete/summary_data/rg_pvals.csv"))
 p.vals.all <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.p.tsv")

 #correlation of factor 1 average
 mean.by.trait <- rowMeans(beta.all[,-1])
 cor.test(mean.by.trait, full.U$U1)
 
#Correlation of V1 with average genetic correlation
 ##Clean up and order the correlation matrix
 name_order <- sapply(ret$trait.names, function(x) which(colnames(rg.mat) == x))
 rg.mat <- rg.mat[name_order, name_order]
 p.mat <- p.mat[name_order, name_order]
 all(colnames(beta) == ret$trait.names)
## Filter to just those signigificant at FDR
 rg.mat.tosum <- rg.mat; diag(rg.mat.tosum) <- 0
 adj.p <- matrix(p.adjust(p=c(p.mat),method="BH"), nrow=137)
 rg.mat.trunc <- rg.mat; diag(rg.mat.trunc) <- 0
 rg.mat.trunc[adj.p > 0.05] <- 0
 ## Make it into a table and joint with V
 sums.by.other <- data.frame("Traits" =colnames(rg.mat.tosum), "rg" = colSums(rg.mat.tosum),
                             "rg_pval_filt"=colSums(rg.mat.trunc))
  rg_and_v <- left_join(full.V, sums.by.other,by="Traits")
  
## Test the correlation 
cor.test(rg_and_v$V1, rg_and_v$rg/137,method ="pearson" )
cor.test(rg_and_v$V1, rg_and_v$rg_pval_filt/137, method="pearson")


##################################
# Number of traits loading on no other factors
which(sapply(apply(ret$V, 1, function(x) which(x != 0)), length) == 1)
##  
  clean.names <- sapply(colnames(p.vals.all)[-1] %>% str_split(., pattern = "\\."), function(x) x[[(length(x)-2)]])
  #Number significant per study withing analysis
  sig.count.df <- data.frame("clean_name" = clean.names, "sig_snps_per" = apply(p.vals.all[,-1], 2, function(x) sum(x < 1e-5, na.rm = TRUE))) %>%
    rename("Traits"=clean_name)
  

##################################  
##Get the ones that are c/c continueous:
  herit.scores <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/panUKBB_analysis/traits_for_review.txt",sep = "#")  %>% separate(sample_size,into = c("case", "control"), sep = "\\/") 
  herit.scores$Neff <- as.numeric(apply(herit.scores, 1, function(x) ifelse(is.na(x[4]), x[3], psych::harmonic.mean(c(as.numeric(x[3]), as.numeric(x[4]))))))
  herit.scores <- herit.scores %>% rowwise() %>% group_by(clean_name) %>% slice_max(Neff) %>% ungroup() %>% arrange(clean_name) %>%
    filter(clean_name %in% ret$trait.names) %>% mutate(clean_name = factor(clean_name, levels = ret$trait.names)) %>% arrange(clean_name) %>%
    distinct() %>% group_by(clean_name) %>% slice_max(Neff, n=1) %>% ungroup() %>% 
    mutate("is_c_c"=ifelse(is.na(control),"cont","cat")) %>% rename("Traits"=clean_name)
  
  
#the number of factors loaded per  trait
  num.loaded.per <- data.frame("Traits" = ret$trait.names,"num_loaded"=sapply(apply(ret$V, 1, function(x) which(x != 0)), length), "V1"=ret$V[,1]) %>% arrange(-1*V1^2) %>%
    left_join(.,herit.scores, by="Traits" )
  
  p.counts.and.loaded <- left_join(num.loaded.per,counts.to.plot %>% rename("Traits"=clean_name),by="Traits" ) %>% left_join(.,sig.count.df,by="Traits")
  top.till.num.loaded<-(p.counts.and.loaded%>% arrange(-V1^2) %>% select(Traits, V1,sig_snps_per,num_loaded,is_c_c))[1:35,]
  17/137
  35/137
  ## Sig SNPs in top V1 traits
  mean(top.till.num.loaded$sig_snps_per)
  ## Sig SNPs in bottom V1 traits
  bottom.num.loaded<-(p.counts.and.loaded%>% arrange(-V1^2) %>% select(Traits, V1,sig_snps_per,num_loaded,is_c_c))[-c(1:35),]
    mean(bottom.num.loaded$sig_snps_per)
  #Go to 35 to have one at 4

#What is the ranking of traits loaded just on V1?
loading.by.rank <- p.counts.and.loaded %>% select(Traits, num_loaded, V1, is_c_c,sig_snps_per, Category.x,Neff.x) %>% mutate("V_score"=V1^2) %>% arrange(-V_score) %>% mutate("rank"=row_number())
    #34 is the top 25
## scree plot has no elbow:
ggplot(loading.by.rank,aes(x=rank,y=abs(V1))) + geom_point()
loading.by.rank %<>% mutate("cos_score"=V1^2/sum(V1^2)) %>% arrange(rank)
loading.by.rank$cum_cos_score <- sapply(1:137, function(i) sum(loading.by.rank$cos_score[1:i]))
#The top 50% of weights
loading.by.rank
sum(filter(loading.by.rank, num_loaded == 1)$cos_score)
loading.by.rank %>% group_by(factor(num_loaded), is_c_c) %>% summarize("count"=n())
### How much do the brain traits contribute?
brain_cat = c("Cognitive", "Neurological", "Psychiatric")
sum((loading.by.rank %>% filter(Category.x %in% brain_cat))$cos_score)#% contribution of these scores
(1/137) * nrow((loading.by.rank %>% filter(Category.x %in% brain_cat))) #average contribution:

ggplot(loading.by.rank, aes(x=cos_score, y=sig_snps_per, color=is_c_c)) +geom_point()
ggplot(loading.by.rank, aes(x=cos_score, y=sig_snps_per, color=Neff.x)) +geom_point()

#Surprisingly, less than the average contribution.
#Get contribution by category
categorical_contributions <- loading.by.rank %>% group_by(factor(Category.x)) %>% summarize("contr_sum"=sum(cos_score), "num_per"=n()) %>%
  mutate("avg_expected"=num_per * (1/137)) %>% mutate("diff_from_avg"=contr_sum - avg_expected)
## Looking at traits loaded just on V1
  #what are the Z-score distributions of these traits that load just on V1, versus others?
  stopifnot(all(se.all$ids == beta.all$ids))
  z.scores.all <- data.frame("ids" = beta.all$ids, as.matrix(beta.all[,-1])/as.matrix(se.all[,-1]))
  colnames(z.scores.all) = c("ids",clean.names)
  z.scores.all %<>% pivot_longer(cols=clean.names,names_to = "Traits") %>% left_join(.,num.loaded.per, by="Traits") %>%
    mutate("v1_only"=ifelse(num_loaded==1, "v1_only","others"))
  ggplot(z.scores.all, aes(x=as.factor(v1_only),y=value^2)) + geom_boxplot()
  summary((z.scores.all %>% filter(num_loaded == 1))$value^2)
  summary((z.scores.all %>% filter(num_loaded > 1))$value^2)
  
  
  #Something better-median per each?
  median.by.cat <- z.scores.all %>% group_by(v1_only,Traits) %>% summarize(median(abs(value),na.rm=TRUE)) %>% set_colnames(c("v1_only", "Trait", "median_Z"))
  ggplot(median.by.cat, aes(x=v1_only, y= median_Z)) + geom_boxplot()
  mean((median.by.cat %>% filter(v1_only == "others"))$median_Z)
  mean((median.by.cat %>% filter(v1_only != "others"))$median_Z)
  
  1-(mean((median.by.cat %>% filter(v1_only != "others"))$median_Z)/mean((median.by.cat %>% filter(v1_only == "others"))$median_Z))  
  
  #Are there lots of brain ones here? No
    #what are the Z-score distributions of these traits that load just on V1, versus others?
  stopifnot(all(se.all$ids == beta.all$ids))
  z.scores.all <- data.frame("ids" = beta.all$ids, as.matrix(beta.all[,-1])/as.matrix(se.all[,-1]))
  colnames(z.scores.all) = c("ids",clean.names)
  z.scores.all %<>% pivot_longer(cols=clean.names,names_to = "Traits") %>% left_join(.,num.loaded.per, by="Traits") %>%
    mutate("v1_only"=ifelse(num_loaded==1, "v1_only","others"))
  ggplot(z.scores.all, aes(x=as.factor(v1_only),y=value^2)) + geom_boxplot()
  summary((z.scores.all %>% filter(num_loaded == 1))$value^2)
  summary((z.scores.all %>% filter(num_loaded > 1))$value^2)
  #Something better-median per each?
  median.by.cat <- z.scores.all %>% group_by(v1_only,Traits) %>% summarize(median(abs(value),na.rm=TRUE)) %>% set_colnames(c("v1_only", "Trait", "median_Z"))
  ggplot(median.by.cat, aes(x=v1_only, y= median_Z)) + geom_boxplot()
  mean((median.by.cat %>% filter(v1_only == "others"))$median_Z)
  mean((median.by.cat %>% filter(v1_only != "others"))$median_Z)
  
  1-(mean((median.by.cat %>% filter(v1_only != "others"))$median_Z)/mean((median.by.cat %>% filter(v1_only == "others"))$median_Z))  

#Are there lots of brain ones here? YES!
num.loaded.per %<>% mutate("brain_related"=ifelse(Category %in% brain_cat, "yes","no"))
#MAKE THE CONTINGENCY TABLE
sing.vs.not <- rbind(table((num.loaded.per %>% filter(num_loaded == 1))$brain_related),table((num.loaded.per %>% filter(num_loaded != 1))$brain_related))
sing.vs.all <- rbind(table((num.loaded.per %>% filter(num_loaded == 1))$brain_related),table((num.loaded.per)$brain_related))
#Is this more than random chance?
fisher.test(t(sing.vs.not))
fisher.test(t(sing.vs.all))
# P-values over the number of traits
## Genomewide significance
sig.per.snp <- apply(p.vals.all[,-1], 1, function(x) sum(x < 5e-8, na.rm=TRUE))
hist(sig.per.snp, xlab="# studies in which a SNP has p<5e-8", main = "Count of traits significantly\nassociated with each SNP")
  #Log version:
  hist.data = hist(sig.per.snp, plot=F)
  hist.data$counts = log(hist.data$counts+1, 10)
  plot(hist.data,xlab="# studies in which a SNP has p<5e-8", main = "Count of traits\nassociated with each SNP", ylab="log10(SNP count)")
by.count <- table(sig.per.snp)

## Nominal significance
sig.per.snp.nominal <- apply(p.vals.all[,-1], 1, function(x) sum(x < 1e-5, na.rm=TRUE))
hist(sig.per.snp.nominal, xlab="# studies in which a SNP has p<1e-5", main = "Count of traits\nassociated with each SNP")
#Log version:
hist.data.n = hist(sig.per.snp.nominal, plot=F)
hist.data.n$counts = log(hist.data.n$counts+1, 10)
plot(hist.data.n,xlab="# studies in which a SNP has p<1e-5", main = "Count of traits\nassociated with each SNP", ylab="log10(SNP count)")
by.count.nominal <- table(sig.per.snp.nominal)
sum(by.count.nominal[10:length(by.count.nominal)])

#save a plot of both
png(file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/suppl_fig_11_snp_sig_count.png")
par(mfrow=c(1,2))
plot(hist.data,xlab="# studies in which\na SNP has p<5e-8", main = "Count of traits\nassociated with each SNP", 
     ylab="log10(SNP count)", yaxp=c(0,5,10))
plot(hist.data.n,xlab="# studies in which a SNP has p<1e-5", main = "Count of traits\nassociated with each SNP", ylab="log10(SNP count)")
dev.off()

##############################
## Factor 1 enrichments (Supp fig 5)
f1.tp <- loadGeneEnrichments(1, "top_fe")
View(f1.tp$enrichments_unique %>% filter(gene_library=="DisGeNET") %>% mutate("custom_by"=p.adjust(`P-value`, method="BH")))

#Factor 1 Tissue enrichments
factor_id <- ret$mapping[1] #get factor 1
tissue.dat <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/NONE_ldsc_enrichment_Multi_tissue_chromatin/plotting_table.csv") %>%
  filter(Source == paste0("F",factor_id)) %>% mutate("Source"="F1")
#MAx tissue across all
tissue.like.eqtl <- tissue.dat %>% group_by(new_category) %>% mutate("category_bfp"=p.adjust(Coefficient_P_value,method = "bonferroni")) %>% slice_min(category_bfp) %>%
  ungroup() %>% mutate("max_tissue_fdr_adjusted"=p.adjust(category_bfp, method="BH"))
#Top Tissue enrichments and Ps 
View((tissue.dat %>% arrange(Factor_specific_fdr))[1:10,])
 

top.10.enrichments <- tissue.dat %>% filter(Source == "F1") %>% arrange(Coefficient_P_value)
color.bars <- top.10.enrichments %>% select(new_category, color_code) %>% distinct()
col.map <- color.bars$color_code; names(col.map) <- color.bars$new_category
ggplot(top.10.enrichments[1:12,], aes(x=reorder(p_tissue,z_score), y=`-log10(P)`, fill=new_category)) +
  geom_bar(stat="identity") + theme_classic(15) + geom_text(aes(label=mark), hjust=0.5) + 
  coord_flip() + xlab("Tissue type") + scale_fill_manual(values=col.map) + labs(fill="Tissue group")+
theme(legend.position = c(0.65, 0.2)) 
ggsave(filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/supp5_enrichment_Bars.png",
       width = 9.5,height=7,units="in")
ggsave(filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/supp5_enrichment_Bars.svg",
       width =9.5,height=7,units="in")

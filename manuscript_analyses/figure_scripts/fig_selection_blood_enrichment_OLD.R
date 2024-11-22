
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")
source("/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/src/get_pLI.R")

pacman::p_load(magrittr, dplyr, ggplot2, data.table, ggsignif)

selection.and.polygen <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/selective_pressure/s_scores.tsv")
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")
#Get blood traits out
blood.trait.idx <- which(ret$V[,4] != 0)
blood.traits <- ret$trait.names[blood.trait.idx]
blood.idx <- which(ret$trait.names %in% blood.traits)
blood.factors <- unique(c(apply(ret$V[blood.idx,], 1, function(x) which(x != 0))))
#Make a nice bar plot
to.bar <- data.frame("traits" = blood.traits, ret$V[blood.trait.idx,blood.factors]) %>% set_colnames(c("Traits", "V1", "V4", "V23", "V11")) %>% pivot_longer(cols=c("V1", "V4", "V23", "V11"))
adj.names <- data.frame("Traits" = unique(to.bar$Traits), "new_names"=c("Mean platelet volume", "Platelet count", "Platelet distr. width"))
to.bar <- left_join(to.bar, adj.names ,by="Traits") %>% mutate("Factor_names"=gsub(name,pattern = "V",replacement="F")) %>%
  mutate("Factor_names"=factor(Factor_names, levels=paste0("F",1:100)))
#Colors:
library(RColorBrewer)
myColors <- brewer.pal(5,"Reds")[c(1,3,5)]
names(myColors) <- unique(to.bar$new_names)
loading.barplot <- ggplot(to.bar %>% filter(name != "V1"), aes(x=value, y=new_names,fill=new_names)) + geom_bar(stat="identity") + facet_wrap(~Factor_names, ncol=1) + theme_classic(16) + geom_vline(xintercept = 0) + 
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position = "right",legend.title =element_blank()) + xlab("V value") + guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values = myColors)

#We want this with a potential enrichment

### see also ~/scratch16-abattle4/ashton/snp_networks/scratch/factor_interpretation/analysis/platelet_development_gene_enrichment.R
#I did it once and lost it. dummy.
## Now for the enrichment test:
#JUST THE IDP GENES
early_Megakaryopoiesis=c("THPO","MPL","GATA1","RUNX1","FLI1","ETV6","GFI1B","HOXA11","MECOM","ANKRD26","RBM8A")
late_mk <- c("AP3B1","HPS1","HPS3","HPS4","HPS5","HPS6","BLOC1S3","BLOC1S6","DTNBP1","LYST","VPS33B","VIPAS39","STXBP2", "NBEA", "NBEAL2", "CYCS","SRC","SLFN14","PLAU","STIM1")
protoplatelet <- c("MYH9", "WAS", "ACTN1", "FLNA", "TUBB1", "DIAPH1", "GP1BA", "GP1BB","GP9","ITGA2B", "ITGB3", "VWF")
platelet_function <- c("P2RY12","TBXA2R", "TBXAS1", "PLA2G4A"," ITGA2B"," RASGRP2","VWF", "ITGB3","FERMTS", "GP1BA","GP9","GP1BB","GP6"," ANO6")

#ADDING IN THE ONES RELATED TO HEREDITARY THROMBOCYTOPENIA (WHERE MISSING)
#OMITTING the grey ones with rare variants causing hereditary HT
em_hmt <- c()
mk_hmt <- c("ABCG5","ABCG8","GNE","MPIG6B")
pp_hmt <- c("ACTB","CDC42")
plt_hmt <- c()

f4.tp <- loadGeneEnrichments(4, "top_fe")
f11.tp <- loadGeneEnrichments(11, "top_fe")
f23.tp <- loadGeneEnrichments(23, "top_fe")


#Permuted and preferred test.
nice.background <- unique(c((f4.tp$all_genes %>% filter(U != 0))$RSID,(f23.tp$all_genes %>% filter(U!=0))$RSID), (f11.tp$all_genes %>% filter(U!=0))$RSID)
full.bg <- (f4.tp$all_genes %>% filter(RSID %in% nice.background) %>% group_by(RSID) %>% slice_head(n=1) %>% ungroup() %>% filter(!is.na(gene)))$gene
#Better names:
factor_sets = list(f4.tp$set_genes, f11.tp$set_genes, f23.tp$set_genes)
factor_titles <- c("F4","F11","F23")
phases=c("Early megakaryopoiesis","Late megakaryopoiesis", "Protoplatelet formation", "Platelet function")
test_sets <- list(early_Megakaryopoiesis,c(late_mk,mk_hmt), c(protoplatelet,pp_hmt),platelet_function)

#Do the enrichment test:
permuted.joint.bg <- tabEnrichmentTestPermuted(factor_sets,full.bg,factor_titles, test_sets,phases,n_permute=10000,conf.level=0.9)
load("/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/enrichment_test_platelets.RData")
#save(permuted.joint.bg,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/enrichment_test_platelets.RData")



factor_sets = list(f4.tp$)
tabEnrichmentTestPermuted(factor_sets,full.bg,factor_titles, sets_to_test,gene_set_names,n_permute=1000,conf.level=0.90)
factor_4_tests <- list(
  testEnrichment((f4.tp$set_genes),full.bg, early_Megakaryopoiesis),
  testEnrichment((f4.tp$set_genes),full.bg, c(late_mk,mk_hmt)),
  testEnrichment((f4.tp$set_genes),full.bg, c(protoplatelet,pp_hmt)), #enriched
  testEnrichment((f4.tp$set_genes),full.bg, platelet_function) )#enriched

factor_23_tests <- list(
  testEnrichment((f23.tp$set_genes),full.bg, early_Megakaryopoiesis),
  testEnrichment((f23.tp$set_genes),full.bg,  c(late_mk,mk_hmt)),
  testEnrichment((f23.tp$set_genes),full.bg, c(protoplatelet,pp_hmt)), #enriched
  testEnrichment((f23.tp$set_genes),full.bg, platelet_function) )#enriched

factor_11_tests <- list(
  testEnrichment((f11.tp$set_genes),full.bg, early_Megakaryopoiesis),
  testEnrichment((f11.tp$set_genes),full.bg, c(late_mk,mk_hmt)),
  testEnrichment((f11.tp$set_genes),full.bg, c(protoplatelet,pp_hmt)), #enriched
  testEnrichment((f11.tp$set_genes),full.bg, platelet_function) )#enriched

joint.background <- rbind(
  data.frame("Factor" = "F23", "p-val" = sapply(factor_23_tests, function(x) x$test$p.value), 
             "OR"=  sapply(factor_23_tests, function(x) x$test$estimate), 
             "LOWER"=sapply(factor_23_tests, function(x) x$test$conf.int[1]),
             "UPPER"=sapply(factor_23_tests, function(x) x$test$conf.int[2]),"Phase"=phases),
  data.frame("Factor" = "F11", "p-val" = sapply(factor_11_tests, function(x) x$test$p.value), 
             "OR"=  sapply(factor_11_tests, function(x) x$test$estimate), 
             "LOWER"=sapply(factor_11_tests, function(x) x$test$conf.int[1]),
             "UPPER"=sapply(factor_11_tests, function(x) x$test$conf.int[2]),"Phase"=phases),
  data.frame("Factor" = "F4", "p-val" = sapply(factor_4_tests, function(x) x$test$p.value), 
             "OR"=  sapply(factor_4_tests, function(x) x$test$estimate), 
             "LOWER"=sapply(factor_4_tests, function(x) x$test$conf.int[1]),
             "UPPER"=sapply(factor_4_tests, function(x) x$test$conf.int[2]),"Phase"=phases)) %>%
  mutate("FDR" = p.adjust(p.val, method="BH")) %>%
  mutate("FDR_sig"=ifelse(FDR < 0.05, "FDR < 0.05", "FDR > 0.05"))
joint.background$Phase <- factor(joint.background$Phase, levels=phases)
joint.background$Factor <- factor(joint.background$Factor, levels = paste0("F",1:100))
ggplot(joint.background, aes(x=Factor, y=OR,color=FDR_sig)) + geom_point(size=3) + 
  geom_errorbar( aes(x=Factor, ymin=LOWER, ymax=UPPER), width=0.3, colour="black", alpha=0.9) + 
  facet_wrap(~Phase, nrow = 1) + theme_classic(16) +
  geom_hline(yintercept = 1,color="black", lty="dashed")+
  scale_color_manual(values=c("red","black")) + theme(legend.title = element_blank()) + ylab("Enrichment OR") + theme(
    legend.position = c(0.1, 0.8), # c(0,0) bottom left, c(1,1) top-right.
    legend.background = element_rect(fill = "white", colour = NA)
  )
#size: 1500 x 400

#Color text by # enerichments
#Plots with just the text for each:
phases=c("Early megakaryopoiesis","Late megakaryopoiesis", "Protoplatelet formation", "Platelet function")
genes.list <- list(early_Megakaryopoiesis,c(late_mk,mk_hmt),c(protoplatelet,pp_hmt),platelet_function)
names(genes.list) <- phases
count.tab <- NULL
for(i in 1:length(genes.list))
{
  g = genes.list[[i]]
  nam = names(genes.list)[i]
  for(gene in g)
  {
    count.tab <- rbind(count.tab, data.frame("set" = nam, "gene"=gene, "in_f11"=sum(f11.tp$set_genes == gene), 
                                             "in_f23"=sum(f23.tp$set_genes == gene),"in_f4"=sum(f4.tp$set_genes == gene)))
  }
  
}

count.tab <- count.tab %>% mutate("col"=case_when(
  in_f11>0 & in_f23>0 & in_f4>0 ~ "gray",
  in_f11>0 & in_f23>0 ~ "green",
  in_f11>0 & in_f4>0 ~ "darkorange",
  in_f4>0 & in_f23>0 ~ "purple",
  in_f11>0 ~ "goldenrod",
  in_f23>0 ~ "cornflowerblue",
  in_f4>0 ~ "red",
  in_f11==0 & in_f23==0 & in_f4==0 ~ "black",
))
count.tab %<>% mutate("font_size"=case_when(
  in_f11>0 & in_f23>0 & in_f4>0 ~ mean(c(in_f11,in_f23,in_f4)),
  in_f11>0 & in_f23>0 ~ mean(c(in_f11,in_f23)),
  in_f11>0 & in_f4>0 ~ mean(c(in_f11,in_f4)),
  in_f4>0 & in_f23>0 ~ mean(c(in_f11,in_f23)),
  in_f11>0 ~ in_f11,
  in_f23>0 ~ in_f23,
  in_f4>0 ~in_f4,
  in_f11==0 & in_f23==0 & in_f4==0 ~ 0,
))

#we want the size to be the # of times its nominated by that
col.list <-unlist(unique(count.tab$col))
names(col.list) <- col.list
# Create the plot
#make x axis:
plt.df <- count.tab %>% filter(font_size > 0) %>% group_by(set) %>% mutate("xcoord"=-1*seq_along(gene)) 
ggplot(plt.df) +
  geom_text(aes(x = 1, y = xcoord, label = gene, color = col, size=font_size+10)) +
  scale_color_manual(values = col.list) +
  theme_void()  + facet_wrap(~set,nrow=1,scales="free_x")



####### Now the blood enrichment:
#Load the gene sets fresh
f23.tp <- loadGeneEnrichments(23, "top_fe")
f11.tp <- loadGeneEnrichments(11, "top_fe")
f4.tp <- loadGeneEnrichments(4, "top_fe")
gsea.f23 <- f23.tp$enrichments %>%  filter(gene_library != "GTEx_Tissue_Expression_Down",gene_library != "GTEx_Tissue_Expression_Up") %>% mutate("FDR"=p.adjust(`P-value`, method="BH"))

gsea.f11 <- f11.tp$enrichments %>%  filter(gene_library != "GTEx_Tissue_Expression_Down",gene_library != "GTEx_Tissue_Expression_Up") %>% mutate("FDR"=p.adjust(`P-value`, method="BH"))

gsea.f4 <- f4.tp$enrichments %>%  filter(gene_library != "GTEx_Tissue_Expression_Down",gene_library != "GTEx_Tissue_Expression_Up") %>% mutate("FDR"=p.adjust(`P-value`, method="BH"))


#keep.groups <- c("GO_Cellular_Component_2023","GO_Molecular_Function_2023", "GO_Biological_Process_2023","KEGG_2021_Human")
keep.groups <- c("GO_Molecular_Function_2023")

#Get just the top 5 per gene group, with FDR < 0.05
top.23.plot <- gsea.f23 %>% filter(FDR < 0.05) %>% arrange(FDR) %>% filter(gene_library %in%  keep.groups) %>%
  group_by(gene_library) %>% slice_head(n=5) %>% ungroup()
top.11.plot <- gsea.f11 %>% filter(FDR < 0.05) %>% arrange(FDR) %>% filter(gene_library %in% keep.groups) %>%
  group_by(gene_library) %>% slice_head(n=5) %>% ungroup()
top.4.plot <- gsea.f4 %>% filter(FDR < 0.05) %>% arrange(FDR) %>% filter(gene_library %in% keep.groups) %>%
  group_by(gene_library) %>% slice_head(n=5) %>% ungroup()


top.in.either <- unique(c(top.4.plot$`Term name`, top.23.plot$`Term name`,top.11.plot$`Term name`))
to.plot <- rbind(gsea.f23 %>% mutate("Source"="F23"), 
                 gsea.f4 %>% mutate("Source"="F4"), gsea.f11 %>% mutate("Source"="F11")) %>% filter(`Term name` %in% top.in.either) %>%
  filter(`Term name` != "Phosphoserine Residue Binding (GO:0050815)") %>% filter(gene_library=="GO_Molecular_Function_2023") %>% mutate("n_genes"=stringr::str_count(Genes,",") + 1)
#Requires being tested in all cases and at least 5 genes contributing to the annotation
keep.terms <- to.plot %>% group_by(`Term name`) %>% mutate("num"=n(), "max_genes"=max(n_genes)) %>% filter(num == 3) %>% filter(max_genes > 5)
#Specify the order
ordered.terms <- c("Transcription Coregulator Binding (GO:0001221)","Transcription Cis-Regulatory Region Binding (GO:0000976)","DNA Binding (GO:0003677)","DNA-binding Transcription Factor Binding (GO:0140297)","Minor Groove Of Adenine-Thymine-Rich DNA Binding (GO:0003680)","GTPase Regulator Activity (GO:0030695)","Tubulin Binding (GO:0015631)","Protein Phosphatase Binding (GO:0019903","Microtubule Binding (GO:0008017)","Cadherin Binding (GO:0045296)")
clean.terms <- gsub(x=ordered.terms, pattern=" \\(GO:\\d+\\)", replacement="")
to.plot$`Term name` <- factor(to.plot$`Term name`, levels = ordered.terms)
to.plot$clean_name <- factor(gsub(x=to.plot$`Term name`, pattern=" \\(GO:\\d+\\)", replacement=""),levels = clean.terms)
enrich.barplot <- ggplot(to.plot %>% filter(`Term name` %in% keep.terms$`Term name`), aes(x=clean_name, y=-log10(`P-value`), fill=Source)) + geom_bar(stat="identity", position="dodge") +
  coord_flip()  + theme_classic(15) + ylab("-log10(p) enrichment") + xlab("GO Molecular function") + theme(legend.position = "bottom")
leg <- ggpubr::get_legend(enrich.barplot)
cowplot::plot_grid(loading.barplot , enrich.barplot,nrow=2,ncol=1)

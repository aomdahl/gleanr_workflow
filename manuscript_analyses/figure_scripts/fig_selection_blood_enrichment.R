LOAD=TRUE
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
clean.blood.names.full <- c("Mean platelet volume", "Platelet count", "Platelet distr. width")
clean.blood.names <- c("MPV", "PCT", "PDW")
to.bar <- data.frame("traits" = blood.traits, ret$V[blood.trait.idx,blood.factors]) %>% set_colnames(c("Traits", "V1", "V4", "V23", "V11")) %>% pivot_longer(cols=c("V1", "V4", "V23", "V11"))
adj.names <- data.frame("Traits" = unique(to.bar$Traits), "new_names"=clean.blood.names, "full_names"=clean.blood.names.full)
to.bar <- left_join(to.bar, adj.names ,by="Traits") %>% mutate("Factor_names"=gsub(name,pattern = "V",replacement="F")) %>%
  mutate("Factor_names"=factor(Factor_names, levels=paste0("F",1:100)))
#Colors:
library(RColorBrewer)
myColors <- brewer.pal(5,"Reds")[c(1,3,5)]
names(myColors) <- unique(to.bar$new_names)
#Colors of factors
factors <- factor(c("F4", "F11", "F23"))
cols <- c("#DC143C", "#DAA520", "#6994F2")
names(cols) <- factors


plot.bars.df <- to.bar %>% filter(name != "V1")
plot.bars.df$Factor_names <- factor(plot.bars.df$Factor_names, levels = c("F4","F11","F23"))
loading.barplot <- ggplot(plot.bars.df, aes(x=value, y=new_names,fill=new_names)) + geom_bar(stat="identity") + facet_wrap(~Factor_names, ncol=1) +
  theme_classic(15) + geom_vline(xintercept = 0) + 
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position = "bottom",legend.title =element_blank(), legend.direction = "horizontal") + 
  xlab("V value") + guides(fill=guide_legend(nrow=1)) + scale_fill_manual(values = myColors)  +
  theme(strip.text = element_text(colour = cols))

ggsave(loading.barplot, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_trait_bars_vert.svg",
       height=7,width=3)
ggsave(loading.barplot, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_trait_bars_vert.png",
       height=7,width=3)
#Opposite orientation:
loading.barplot.long <- ggplot(to.bar %>% filter(name != "V1"), aes(x=value, y=new_names,fill=new_names)) + geom_bar(stat="identity") + facet_wrap(~Factor_names, nrow=1) +
  theme_classic(15) + geom_vline(xintercept = 0) + 
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position = "bottom",legend.title =element_blank(), legend.direction = "horizontal") + xlab("V value") +
  scale_fill_manual(values = myColors)
ggsave(loading.barplot.long, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_trait_bars.svg",
       height=7,width=3)
ggsave(loading.barplot.long, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_trait_bars.png",
       height=7,width=3)


#FOR POSTER VERSINO:
ggplot(to.bar %>% filter(name != "V1"), aes(x=value, y=full_names,fill=new_names)) + geom_bar(stat="identity") + facet_wrap(~Factor_names, nrow=1) +
  theme_classic(17) + geom_vline(xintercept = 0) + 
  theme(legend.position = "none",legend.title =element_blank(), legend.direction = "horizontal") + xlab("V value") +
  scale_fill_manual(values = myColors) + ylab("")
ggsave(filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_trait_bars.POSTER.svg",
       height=2,width=9)
ggsave( filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_trait_bars.POSTER..png",
       height=2,width=9)
### see also ~/scratch16-abattle4/ashton/snp_networks/scratch/factor_interpretation/analysis/platelet_development_gene_enrichment.R


#I did it once and lost it. dummy.
## Now for the enrichment test:
#JUST THE IDP GENES
early_Megakaryopoiesis=c("THPO","MPL","GATA1","RUNX1","FLI1","ETV6","GFI1B","HOXA11","MECOM","ANKRD26","RBM8A")
late_mk <- c("AP3B1","HPS1","HPS3","HPS4","HPS5","HPS6","BLOC1S3","BLOC1S6","DTNBP1","LYST","VPS33B","VIPAS39","STXBP2", "NBEA", "NBEAL2", "CYCS","SRC","SLFN14","PLAU","STIM1")
protoplatelet <- c("MYH9", "WAS", "ACTN1", "FLNA", "TUBB1", "DIAPH1", "GP1BA", "GP1BB","GP9","ITGA2B", "ITGB3", "VWF")
platelet_function <- c("P2RY12","TBXA2R", "TBXAS1", "PLA2G4A"," ITGA2B"," RASGRP2","VWF", "ITGB3","FERMTS", "GP1BA","GP9","GP1BB","GP6"," ANO6")
macro.genes <- c("THPO","ANKRD26","ETV6","FLI1","GATA1","GFI1B","HOXA11","MECOM","RUNX1","NBEAL2","ABCG5","ABCG8","GNE","MPIG6B","SLFN14","SRC","ACTB","ACTN1","CDC42","DIAPH1","FLNA","MYH9","TUBB1","GP1BA","GP1BB","GP9","ITGA2B","ITGB3","VWF")

#ADDING IN THE ONES RELATED TO HEREDITARY THROMBOCYTOPENIA (WHERE MISSING)
#OMITTING the grey ones with rare variants causing hereditary HT
em_hmt <- c()
mk_hmt <- c("ABCG5","ABCG8","GNE","MPIG6B")
pp_hmt <- c("ACTB","CDC42")
plt_hmt <- c()

f4.tp <- loadGeneEnrichments(4, "top_fe")
f11.tp <- loadGeneEnrichments(11, "top_fe")
f23.tp <- loadGeneEnrichments(23, "top_fe")

#For viz- the scree plots of each?



#Permuted and preferred test.
nice.background <- unique(c((f4.tp$all_genes %>% filter(U != 0))$RSID,(f23.tp$all_genes %>% filter(U!=0))$RSID), (f11.tp$all_genes %>% filter(U!=0))$RSID)
full.bg <- (f4.tp$all_genes %>% filter(RSID %in% nice.background) %>% group_by(RSID) %>% slice_head(n=1) %>% ungroup() %>% filter(!is.na(gene)))$gene
#Better names:
factor_sets = list(f4.tp$set_genes, f11.tp$set_genes, f23.tp$set_genes)
factor_titles <- c("F4","F11","F23")
or_test_hmt <- testEnrichment(unique(f4.tp$set_genes), unique(full.bg), macro.genes,conf.level=0.90,alternative="two.sided")
plot.me <- data.frame("OR"=or_test_hmt$test$estimate, "upper"=or_test_hmt$test$conf.int[2], "lower"=or_test_hmt$test$conf.int[1], "p"=or_test_hmt$test$p.value, "Factor"="U4")
#Quick nice plot of this 10/28
ggplot(plot.me, aes(x=Factor, y=OR)) + geom_point(size=3) + 
  geom_errorbar( aes(x=Factor, ymin=lower, ymax=upper), width=0.3, colour="black", alpha=0.9) + theme_classic(16) +
  geom_hline(yintercept = 1,color="black", lty="dashed")+ ylab("Enrichment OR")  + coord_flip()




##Harder test- just within blood genes
blood.bg <- c(f4.tp$set_genes,f11.tp$set_genes,f23.tp$set_genes)
testEnrichment(unique(f4.tp$set_genes), unique(blood.bg), macro.genes,conf.level=0.90,alternative="two.sided")
testEnrichment(unique(f11.tp$set_genes), unique(blood.bg), macro.genes,conf.level=0.90,alternative="two.sided")
testEnrichment(unique(f23.tp$set_genes), unique(blood.bg), macro.genes,conf.level=0.90,alternative="two.sided")






## Next test- all IDP genes at stages
phases=c("Early megakaryopoiesis","Late megakaryopoiesis", "Protoplatelet formation", "Platelet function")
test_sets <- list(early_Megakaryopoiesis,c(late_mk,mk_hmt), c(protoplatelet,pp_hmt),platelet_function)
#tabEnrichmentTestPermuted <- function(factor_sets,full.bg,factor_titles, sets_to_test,gene_set_names,n_permute=1000)
if(LOAD)
{
  load("/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/enrichment_test_platelets.RData")
  
}else
{
  ##First test- just those for macro genes
  permuted.joint.hmtc <- tabEnrichmentTestPermuted(factor_sets,full.bg,factor_titles, list(macro.genes),c("HMTC"),n_permute=10000,conf.level=0.9)
  permuted.joint.bg <- tabEnrichmentTestPermuted(factor_sets,full.bg,factor_titles, test_sets,phases,n_permute=50000,conf.level=0.9)
  save(permuted.joint.bg,permuted.joint.hmtc,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/enrichment_test_platelets.RData")
  
}


permuted.joint.bg$Phase <- factor(permuted.joint.bg$Phase, levels=phases)
permuted.joint.bg$Factor <- factor(permuted.joint.bg$Factor, levels = paste0("F",1:100))
permuted.joint.bg$fdr_sig <- ifelse(permuted.joint.bg$FDR < 0.05, "FDR < 0.05", "FDR > 0.05")



#size: 1500 x 400
custom.markers <- ggplot(permuted.joint.bg, aes(x=Factor, y=OR,color=-log10(p.val))) + geom_point(size=3) + 
  geom_errorbar( aes(x=Factor, ymin=LOWER, ymax=UPPER), width=0.3, colour="black", alpha=0.9) + 
  facet_wrap(~Phase, nrow = 1) + theme_classic(18) +
  geom_hline(yintercept = 1,color="black", lty="dashed") + scale_color_gradient(low="blue",high="red",limits=c(0,5)) + 
  labs(color="Permutation\n-log10(q)") + ylab("Enrichment OR") + theme(
    legend.position = "none", # c(1, 0.8) bottom left, c(1,1) top-right.
    legend.background = element_rect(fill = "transparent", colour = NA),
                                     legend.text = element_text(size=12),
                                     legend.title=element_text(size=12), 
                                     legend.title.position = "right",  plot.title =element_text(size=10, face='bold', color="grey"),
    axis.title.x = element_blank(),axis.text.x= element_text(colour=cols, size=18)) #+
  #ggtitle("IPD genes by differention stage")
ggsave(custom.markers, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_idp_enrichment.svg",
        height=4,width=14)

ggsave(custom.markers, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_idp_enrichment.png",
       height=4,width=15)
#size: 1500 x 400



####### Now the blood enrichment:


###Panglao DB:
all.markers <- fread("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")
#unique(all.markers$`cell type`)
platelet.genes <- (dplyr::filter(all.markers, `cell type` == "Platelets"))$`official gene symbol`
hsc.genes <- (dplyr::filter(all.markers, `cell type` == "Hematopoietic stem cells"))$`official gene symbol`
mk.genes <- (dplyr::filter(all.markers, `cell type` == "Megakaryocytes"))$`official gene symbol`

factor_sets = list(f4.tp$set_genes, f11.tp$set_genes, f23.tp$set_genes)
factor_titles <- c("F4","F11","F23")
sets_to_test <- list(platelet.genes,hsc.genes,mk.genes)
gene_set_names <- c("Platelets","HSCs", "Megakaryocytes")
if(LOAD)
{
  load("/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/enrichment_test_panglao.RData")
}else
{
  unique.permuted.cell.types <- tabEnrichmentTestPermuted(factor_sets,f11.tp$bg_genes,factor_titles, sets_to_test,gene_set_names,n_permute=10000,conf.level=0.9)
  save(unique.permuted.cell.types,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/enrichment_test_panglao.RData")
}

unique.permuted.cell.types$Phase <- factor(unique.permuted.cell.types$Phase, levels=c("HSCs", "Megakaryocytes","Protoplatelets", "Platelets"))
unique.permuted.cell.types$Factor <- factor(unique.permuted.cell.types$Factor, levels = paste0("F",1:100))

panglao.markers <- ggplot(unique.permuted.cell.types, aes(x=Factor, y=OR,color=-log10(p.val))) + geom_point(size=3) + 
  geom_errorbar( aes(x=Factor, ymin=LOWER, ymax=UPPER), width=0.3, colour="black", alpha=0.9) + 
  facet_wrap(~Phase, nrow = 1,drop=FALSE) + theme_classic(18) +
  geom_hline(yintercept = 1,color="black", lty="dashed")+ ylab("Enrichment OR")  + 
  labs(color="Permutation\n-log10(p)") + scale_color_gradient(low="blue",high="red",limits=c(0,5))  + 
  theme(
    legend.position = "left", # c(1, 0.8) bottom left, c(1,1) top-right.
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(size=12),
    legend.title=element_text(size=12), 
    legend.title.position = "right", plot.title =element_text(size=10, face='bold', color="grey"),
    axis.text.x= element_text(colour=cols, size=18))# + ggtitle("Cell type gene marker enrichment")
ggsave(panglao.markers, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_panglao_enrichment.svg",
       height=4,width=15)


#Move the legend to the empty spot
fitted_plot <- ggplot(unique.permuted.cell.types, aes(x=Factor, y=OR,color=-log10(p.val))) + geom_point(size=3) + 
  geom_errorbar( aes(x=Factor, ymin=LOWER, ymax=UPPER), width=0.3, colour="black", alpha=0.9) + 
  facet_wrap(~Phase, nrow = 1,drop=FALSE) + theme_classic(18) +
  geom_hline(yintercept = 1,color="black", lty="dashed")+ ylab("Enrichment OR")  + 
  labs(color="Permutation\n-log10(p)") + scale_color_gradient(low="blue",high="red",limits=c(0,5))  + 
  theme(
    legend.position =c(0.65, 0.6), #bottom left, c(1,1) top-right.
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(size=12),
    legend.title=element_text(size=12), 
    legend.title.position = "right", plot.title =element_text(size=10, face='bold', color="grey"),
    axis.text.x= element_text(colour=cols,size=18)) #+ ggtitle("Cell type gene marker enrichment")

panglao.fitted <- remove_facets(fitted_plot, "aa#a")

#size: 1500 x 400

cowplot::plot_grid(custom.markers, panglao.markers,ncol=1)
ggsave(filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_both_enrichments.svg",
       height=5,width=13)

#This is the one I want.
cowplot::plot_grid(custom.markers , panglao.fitted,ncol=1, rel_heights = c(0.9,1))
ggsave(filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig4_both_enrichments_with_blank.svg",
       height=7,width=13)



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


#for the paper-specific results, see `blood_factor_enrichr_results.Rmd`
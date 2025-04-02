### Fig4 description of S
pacman::p_load(magrittr, dplyr, ggplot2,data.table)
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")

##Factor 1 sparsity
cat(round(100*sum(ret$U[,1]==0)/nrow(ret$U), digits=2), "% of SNPs in factor 1 have a value of 0\n")
##Factors 32 and 18 are weakly correlated:
cat("Factors 32 and 18 have a correlation of", cor.test(ret$U[,32], ret$U[,18])$estimate, "\n")

#Selection scores:
sel.dat <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final//selective_pressure/s_scores.tsv")
cat("Selection score of factor 1 is:", round((sel.dat %>% filter(Factor=="U1"))$S_hat,digits=2), "\n")
cat("Selection score of factor 32 is:", round((sel.dat %>% filter(Factor=="U32"))$S_hat,digits=2), "\n")
cat("Minimum selection score is:", round(min((sel.dat)$S_hat),digits=2)," in factor",(sel.dat %>% arrange(S_hat))[1,]$Factor,"\n")

#Enrichments in tissue markers
#Note- have to use mapping here because factors were out of ordered on run (bug in PVE code.)
#This has ince been fixed, with a mapping between the new and the old factors.
tissue.plotting.dat <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/NONE_ldsc_enrichment_Multi_tissue_chromatin//plotting_table.csv")
mapping <- data.frame("new" = 1:58, "old" = ret$mapping)
mapping.df <- mapping %>% mutate("Source"=paste0("F",old)) %>% mutate("Source_new"=paste0("F",new))
#For refrence of all
labels <- labelsPrep()
P="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/NONE_ldsc_enrichment_Multi_tissue_chromatin/"
f18.all.tissue <- fread(paste0(P, (mapping.df %>% filter(Source_new=="F18"))$Source,".multi_tissue.cell_type_results.txt")) %>% left_join(.,labels, by="Name")
f32.all.tissue <- fread(paste0(P, (mapping.df %>% filter(Source_new=="F32"))$Source,".multi_tissue.cell_type_results.txt")) %>% left_join(.,labels, by="Name")#F28

#Get top tissue type enrichments by category
#MAx tissue across all
# <- evalTissuesLikeEQTLs(tissue.plotting.dat, "F18") %>% arrange(max_tissue_fdr_adjusted) %>% print()
#f32.like.eqtl <- evalTissuesLikeEQTLs(tissue.plotting.dat, "F32") %>% arrange(max_tissue_fdr_adjusted) %>% print()
#A better version of this:
f18.like.eqtl <- evalTissuesLikeEQTLs(f18.all.tissue) %>% arrange(max_tissue_fdr_adjusted) %>% print()
f32.like.eqtl <- evalTissuesLikeEQTLs(f32.all.tissue) %>% arrange(max_tissue_fdr_adjusted) %>% print()


  tissue.plotting.dat <- left_join(tissue.plotting.dat, mapping.df, by="Source")
  tissue.plotting.dat$Source <- tissue.plotting.dat$Source_new
  (tissue.plotting.dat %>% filter(Source == "F32") %>% arrange(factor_tissue_fdr) %>% 
    select(Name, Source, tissue,category, global_fdr, p_tissue, Factor_specific_fdr))[1:5,]
  (tissue.plotting.dat %>% filter(Source == "F18") %>% arrange(factor_tissue_fdr)%>% 
    select(Name, Source, tissue,category, global_fdr, p_tissue, Factor_specific_fdr))[1:5,]


#Look at gene enrichments
f32.tp <- loadGeneEnrichments(32, "top_fe")
(f32.tp$enrichments_unique %>% filter(gene_library == "Human_Phenotype_Ontology"))[1:10,]
f32.tp$enrichments_unique %>% filter(gene_library == "DisGeNET")
f32.tp$enrichments_unique %>% filter(gene_library == "KEGG_2021_Human")
f32.tp$enrichments_unique %>% filter(gene_library == "Chromosome_Location")

f18.tp <- loadGeneEnrichments(18, "top_fe")
(f18.tp$enrichments_unique %>% filter(gene_library == "Human_Phenotype_Ontology"))[1:10,]
f18.tp$enrichments_unique %>% filter(gene_library == "DisGeNET")

#Look at SNPs:
## Load SNP data
full.snps <- loadMAFData()

full.U <- data.frame("SNPs"=ret$snp.ids, ret$U) %>% magrittr::set_colnames(c("SNPs", paste0("U", 1:ret$K))) 
extreme_scores <- full.U %>% select(SNPs, U32, U1, U18, U37) %>% left_join(., full.snps %>% rename("SNPs"=rsid), by = "SNPs") 


### SNP selection plot
for.maf.plot <- extreme_scores %>% mutate("F1"=scale(U1,center = FALSE), "F32"=scale(U32, center=FALSE)) %>%
  select(SNPs, F1,F32, maf) %>%
  pivot_longer(cols=c(F1,F32),names_to = "source_col",values_to = "U_scores") %>% 
  mutate("source_col"=ifelse(source_col == "F1", "F1 (ubiquitous)", source_col))  %>%
  set_colnames(c("SNPs", "maf","source_col","U_scores"))

dot.maf <- ggplot(for.maf.plot, aes(y=U_scores^2, x=maf, color = maf)) + geom_point() +
  scale_color_gradient(low="yellow",high="black") + facet_wrap(~source_col) + theme_bw(14) + 
  xlab("SNP minor allele frequency") + ylab(bquote("Scaled SNP loading (" *U^2 *")")) + 
                                              labs(color="MAF") +  theme(strip.background =element_rect(fill="white"))
                                            #density of the points is misleading
ggsave(dot.maf,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_maf_dotplot.svg", width=6, height=2)
ggsave(dot.maf,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_maf_dotplot.png", width=6, height=2)

#Dot maf vertical:
#dot.maf.v <- ggplot(for.maf.plot, aes(y=U_scores^2, x=maf, color = maf)) + geom_point() +
#  scale_color_gradient(low="skyblue",high="black") + facet_wrap(~source_col,ncol=1) + theme_bw(14) + 
#  xlab("Minor allele frequency") + ylab(bquote("Scaled SNP loading (" *U[k]^2 *")")) + 
#  labs(color="MAF") +  theme(strip.background =element_rect(fill="white"), legend.position = "bottom", 
#                             legend.key.height = unit(0.5, 'cm'),
#                             legend.key.width = unit(1.7,"cm"))
dot.maf.v <- ggplot(for.maf.plot, aes(y=U_scores^2, x=maf)) + geom_point() + facet_wrap(~source_col,ncol=1) + theme_bw(15) + 
  xlab("Minor allele frequency") + ylab(bquote("Scaled SNP loading (" *U[k]^2 *")")) + 
  labs(color="MAF") +  theme(strip.background =element_rect(fill="white"))
#density of the points may be confusing
ggsave(dot.maf.v,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_maf_dotplot_vertical.svg", width=4, height=6)
ggsave(dot.maf.v,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_maf_dotplot_vertical.png", width=4, height=6)
#For poster version
ggsave(dot.maf.v,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_maf_dotplot_vertical.POSTER.svg", width=4.25, height=4.41)
ggsave(dot.maf.v,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_maf_dotplot_vertical.POSTER.png", width=4.25, height=4.41)
#Density plot....

hist.poly <- ggplot(for.maf.plot %>% filter(U_scores != 0), aes(x=U_scores,fill=source_col)) + 
  geom_histogram(position="identity", alpha=0.5,bins=200) + theme_classic(15)+
  xlim(c(-5,5)) + xlab("Non-zero SNP loadings (U)") + ylab("# SNPs") + 
  labs(fill="Factor") + theme(legend.position = "inside",
    legend.position.inside = c(0.8,0.5 )) 

#I think we are dropping the hist
cowplot::plot_grid(dot.maf, hist.poly)
ggsave(file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_selection_maf.svg")

#Mini image for Figure 1 thumbnail:
ggplot(for.maf.plot %>% filter(source_col == "F32"), aes(y=U_scores^2, x=maf)) + geom_point() + theme_classic(15) + 
theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
ggsave(file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_1_selection_maf.svg", height = 4,width=5)
ggsave(file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_1_selection_maf.png", height = 4,width=5)
#######################the nice dotplot##
## Load in the data
selection.and.polygen <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/selective_pressure/s_scores.tsv")
num_loaded.by.factor <- data.frame("Factors" = paste0("V", 1:ret$K), "num_loaded" = apply(ret$V,2, function(x) sum(x != 0))) %>% mutate("designation"=ifelse(num_loaded == 1, "singleton", ifelse(num_loaded == 137, "ubiquitous", "std")))

count_and_polygen <- left_join(selection.and.polygen, num_loaded.by.factor %>% mutate("Factor"=gsub(x=Factors, pattern="V", replacement="U")), by="Factor")
count_and_polygen %<>% mutate("num_factors_f" = ifelse(num_loaded > 100, "15+", as.character(num_loaded))) %>%
  mutate("num_factors_f"=factor(num_factors_f, levels =c("1","2","3","4","5","6","7","13","15+"))) %>%
  mutate("num_factors_c" = ifelse(num_loaded > 100, 18, (num_loaded))) 
##THIS ONE
library(RColorBrewer)
#old version- discrete scale
#dotplot.top <- ggplot(count_and_polygen , aes(x=Me, y=S_hat,fill = num_factors_f)) +
#  geom_point(size=3, shape=21, color="black", stroke=0.5) + xlab(bquote("Polygenicity ("~ M[e] ~")")) + 
#  ylab(bquote(hat(S))) + theme_bw(14) + labs(fill="# of\nloaded\ntraits") + scale_fill_brewer(palette = "Blues") +
#  theme(legend.position = "inside", legend.position.inside = c(0.6,0.2),
#        legend.direction="horizontal",legend.box.background = element_rect(colour = "black"))
dotplot.top <-  ggplot(count_and_polygen , aes(x=Me_scaled, y=S_hat,fill = num_factors_c)) +
  geom_point(size=3, shape=21, color="black", stroke=0.5) + xlab(bquote("Polygenicity")) + 
  ylab(bquote("Selection signature ("*hat(S)*")")) + theme_bw(14) + labs(fill="# of\nloaded\ntraits") + 
  scale_fill_gradient(low=brewer.pal(9,"Reds")[1], high=brewer.pal(9,"Reds")[8],
                      breaks = c(2,6,10,14,18),labels = c("2","6","10","14","18+")) +
  theme(legend.position = "inside", legend.position.inside = c(0.6,0.2),
        legend.direction="horizontal",legend.box.background = element_rect(colour = "black"),
        legend.key.height = unit(0.5, 'cm'),legend.key.width = unit(1, 'cm'), legend.margin = margin(5, 14, 5, 13))

#For poster: 723x 340
# labels = c("low", "med", "high")
ggsave(dotplot.top,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_selection_s_poly_dots.svg",width = 7,height=5)
ggsave(dotplot.top,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_selection_s_poly_dots.png",width = 7,height=5)
#For the poster
ggsave(dotplot.top,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_selection_s_poly_dots.POSTER.svg",width = 8.5,height=4)
ggsave(dotplot.top,file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_selection_s_poly_dots.POSTER.png",width = 8.5,height=4)


cowplot::plot_grid(dotplot.top, dot.maf,labels=c("A","B"))
#try 1400 gby 400
ggsave(file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_selection_dot_maf.svg",width = 13,height=2.5)
ggsave(file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_selection_dot_maf.png",width = 13,height=2.5)

cowplot::plot_grid(dotplot.top,dot.maf.v,rel_widths = c(1,0.6))
ggsave(file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_selection_dot_maf_vertical.svg",width = 10,height=4)
ggsave(file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig_4_selection_dot_maf_vertical.png",width = 10,height=4)

#Color the spots to manually identify in inkscape:
col.dots <- c("U32", "U37", "U43", "U18", "U11", "U4", "U23", "U1")
count_and_polygen.color <- count_and_polygen %>% mutate("fill"=ifelse(Factor %in% col.dots, "red", "black"))
ggplot(count_and_polygen.color , aes(x=Me_scaled, y=S_hat,fill = fill)) +
  geom_point(size=3, shape=21, color="black", stroke=0.5) + xlab(bquote("Polygenicity")) + 
  ylab(bquote("Selection signature ("*hat(S)*")")) 

#Color by the average sample size to ensure our estimates aren't biased that way.
ret$trait.names <- ret$trait.names %>% gsub(x=., pattern = "X_", replacement = "_") %>% gsub(x=., pattern = "non.albumin_protein", replacement = "non-albumin_protein")
herit.scores <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/panUKBB_analysis/traits_for_review.txt",sep = "#")  %>% separate(sample_size,into = c("case", "control"), sep = "\\/") 
herit.scores$Neff <- as.numeric(apply(herit.scores, 1, function(x) ifelse(is.na(x[4]), x[3], psych::harmonic.mean(c(as.numeric(x[3]), as.numeric(x[4]))))))
herit.scores %<>% rowwise() %>% group_by(clean_name) %>% slice_max(Neff) %>% ungroup() %>% arrange(clean_name) %>%
  filter(clean_name %in% ret$trait.names) %>% mutate(clean_name = factor(clean_name, levels = ret$trait.names)) %>% arrange(clean_name) %>%
  distinct() %>% group_by(clean_name) %>% slice_max(Neff, n=1) %>% ungroup() %>% 
  mutate("is_c_c"=ifelse(is.na(control),"cont","cat")) %>% rename("Traits"=clean_name)
with.trait.dat <- data.frame("Traits" = ret$trait.names, apply(ret$V, 2, function(x) x^2/sum(x^2))) %>% 
  left_join(.,herit.scores %>% select(Traits, Neff) %>% distinct(), by = "Traits") %>% distinct()
colnames(with.trait.dat)
per.factor.avg.neff <- apply(with.trait.dat[,c(2:59)],2, function(x) sum(x*with.trait.dat$Neff))
count_and_polygen.color$factor_n <- per.factor.avg.neff
par(mfrow=c(1,2))
plot(count_and_polygen.color$factor_n,count_and_polygen$Me_scaled, col="red", pch=19, xlab='Factor N', ylab="Polygenicity")
plot(count_and_polygen.color$factor_n,count_and_polygen$S_hat, col="blue", pch=19, xlab='Factor N', ylab="S")
cor.test(count_and_polygen.color$factor_n,count_and_polygen$Me_scaled)
cor.test(count_and_polygen.color$factor_n,count_and_polygen$Me_scaled)
#########################
#Scree plots of factor distributions
scree.f4.df <- data.frame("U4_sq" = ret$U[,4]^2) %>% arrange(-U4_sq) %>% mutate("rank"=row_number())
ggplot(scree.f4.df,aes(x=rank,y=U4_sq)) + geom_jitter() + theme_classic()
y_range <- range(my_data$y)
p2_dens <- ggplot(scree.f4.df, aes(x = U4_sq)) +
  geom_density() +scale_x_continuous(trans = 'log10',limits =c(0,6.5e-6 )) +coord_flip()

  sorted.valsdf <- data.frame("x"=rank(-ret$U[,4]^2), "y"=ret$U[,4]^2) %>% arrange(x)
  elbow.i <- ceiling(pathviewr::find_curve_elbow(sorted.valsdf, plot_curve = TRUE) * (1/3))
  ggplot(sorted.valsdf,aes(x=x,y=y)) + geom_jitter() + theme_classic(16)+ 
    geom_hline(yintercept =  sorted.valsdf[elbow.i, 2 ], lty="dashed",col="blue") + 
    xlab(bquote("SNP" ~ U[4] ~"rank")) + ylab(bquote(U[4]^2))
nrow(ret$U ) - elbow.i

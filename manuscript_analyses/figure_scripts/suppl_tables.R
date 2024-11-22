############################################################################################################
# Ashton Omdahl, November 2024
### Supplementary table generation:
############################################################################################################
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")
source("/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/src/get_pLI.R")

pacman::p_load(magrittr, dplyr, ggplot2, data.table, ggsignif,RColorBrewer)
# Supplementary table 1- phenotypes included in UKBB/Finngen analysis:

#finngen download: /data/abattle4/lab_data/GWAS_summary_statistics/FinnGen/finngen_download_commands_jan2024.sh
# Supplementary table  2:
load("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")
## Get groups
group.assigns <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/panUKBB_analysis/traits_for_review.txt", sep="#")
## Get trait names
traits.to.verify <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete.trait_list.tsv",header = FALSE)
## Get the clean trait names
clean_trait_names <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/trait_clean_names.csv",header = FALSE) %>%
  select(V1,V2) %>% set_colnames(c("clean_name", "cleaner"))

counts.to.plot <- filter(group.assigns, clean_name %in% traits.to.verify$V1) %>% group_by(clean_name) %>% slice_head(n=1) %>%
  ungroup() %>% left_join(., clean_trait_names, by="clean_name")


herit.scores <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/panUKBB_analysis/traits_for_review.txt",sep = "#")  %>% separate(sample_size,into = c("case", "control"), sep = "\\/") 
herit.scores$Neff <- as.numeric(apply(herit.scores, 1, function(x) ifelse(is.na(x[4]), x[3], psych::harmonic.mean(c(as.numeric(x[3]), as.numeric(x[4]))))))
herit.scores <- herit.scores %>% rowwise() %>% group_by(clean_name) %>% slice_max(Neff) %>% ungroup() %>% arrange(clean_name) %>%
  filter(clean_name %in% ret$trait.names) %>% mutate(clean_name = factor(clean_name, levels = ret$trait.names)) %>% arrange(clean_name) %>%
  distinct() %>% group_by(clean_name) %>% slice_max(Neff, n=1) %>% ungroup() %>% 
  mutate("is_c_c"=ifelse(is.na(control),"cont","cat")) 
clean.combined <- left_join(counts.to.plot,herit.scores, by= "clean_name") %>% select(cleaner,clean_name, description.x,Category.x, case, control,Neff.y, is_c_c, estimates.final.h2_liability.x) %>%
  set_colnames(c("Trait","clean_name", "Description", "Category", "Num cases", "Num controls", "Neff", "Continuous_or_categorical", "PanUKBB_h2_liability_estimate")) %>% 
  select(-Continuous_or_categorical)
#add the download call?
trait.cals <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/panUKBB_complete.studies.tsv")
lookups <- basename(trait.cals$V1) %>% gsub(x=., pattern=".EUR_ONLY.gz", replacement = "")
phenotype.names <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/Pan-UK_Biobank_phenotype_manifest-phenotype_manifest.tsv", sep = "\t") %>%
  filter(any(sapply(lookups, function(y) grepl(pattern = y, x=wget))))

keep.calls <- sapply(lookups, function(l) which(grepl(pattern = l, x=phenotype.names$wget)))
panukk.meta.dat <- cbind(trait.cals, phenotype.names[keep.calls,c(1:7,79)]) %>% 
  set_colnames(c("Local_path", "clean_name", "genome_build", "Neff", "Cohort", "Trait_type","phenocode", "Sample_sex", "Coding","Coding_modifier", "Description","Description_more", "wget")) %>%
  select(-Neff, -Local_path, -genome_build, -Cohort, -Description_more, -Description)

final.list <- left_join(clean.combined, panukk.meta.dat,by="clean_name") %>% mutate("Coding" = paste0(Coding, "_", Coding_modifier)) %>%
  select(Trait, clean_name, Trait_type, phenocode, Coding, Sample_sex, Category, `Num cases`, `Num controls`, Neff, Description, PanUKBB_h2_liability_estimate,wget )
write.table(final.list, file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/Suppl_tab_2.tsv",quote = FALSE, row.names = FALSE,sep="\t")

# Supplementary table 3: All tissue enrichments data.

tissue.dat <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/NONE_ldsc_enrichment_Multi_tissue_chromatin/plotting_table.csv")  %>%
  select(Name, tissue, mark, new_category,Source,Coefficient, Coefficient_std_error, Coefficient_P_value, global_fdr, Factor_specific_fdr) %>%
  rename("Category"=new_category, "Factor"=Source)
mapping <- data.frame("new" = 1:58, "old" = ret$mapping)
mapping.df <- mapping %>% mutate("Factor"=paste0("F",old)) %>% mutate("Factor_new"=paste0("F",new))
tissue.dat <- left_join(tissue.dat, mapping.df, by="Factor")
tissue.dat$Factor <- tissue.dat$Factor_new

write.table(tissue.dat %>% select(-new,-old,-Factor_new) %>% arrange(Factor), file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/Suppl_tab_3.csv",quote = FALSE, row.names = FALSE,sep=",")
##############
# Supplementary Table 15: Simulation data:
#Save the data with different combinations of overlap for each trait
#the overlap list used in our simulations:
overlap.list <- c(0,10000,15000,20000,23000,26000,30000)
all.options <- expand.grid(overlap.list, overlap.list, overlap.list)
all.ids <- sapply(1:nrow(all.options), function(i)  paste(sort(unlist(all.options[i,])), collapse="-"))
sub.options <- all.options[!duplicated(all.ids),] %>% mutate("total"=Var1+Var2+Var3) %>% mutate("tot_prop"=total/90000)
iter=10
count.per.overlap.level <- sub.options %>% group_by(prop) %>% summarize("num_sims"=n())
ggplot(count.per.overlap.level, aes(x=prop,y=num_sims*iter)) + geom_bar(stat="identity") + theme_classic() + 
  xlab("Proportion of all 90,000 simulated samples shared") +
  ylab("Number of simulations performed")

n=30000
write.sim.options <- sub.options %>% mutate("Prop_cohort1"=Var1/n, "Prop_cohort2"=Var2/n, "Prop_cohort3"=Var3/n) %>%
  select(Prop_pop1,Prop_pop2,Prop_pop3, total ) %>% mutate("Total_shared_prop"=total/(n*3)) %>%
  select(-total) %>% arrange(Total_shared_prop)
write.table(write.sim.options, file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/Suppl_tab_15.csv",quote = FALSE, row.names = FALSE,sep=",")


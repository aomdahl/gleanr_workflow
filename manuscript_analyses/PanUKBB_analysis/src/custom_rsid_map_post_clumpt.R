#!/usr/bin/env Rscript
pacman::p_load(magrittr, dplyr, data.table)
args = commandArgs(trailingOnly=TRUE)
#args[1] <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.250kb.0.2r2.prune.in"
#args[2] <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.pruned_rsids.txt"
setwd("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization")
ref.lookup <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz") 
#ids.reference <- fread("gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.pruned_rsids.intermediate.txt",header = FALSE)
#prune.set <- fread("gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.250kb.0.2r2.prune.in")
prune.set <- fread(args[1])
#Maybe some are redundant?
#Just do it straight, none of this generalizeable garbage. No more time for that ish.
sub.down <- ref.lookup %>% mutate("snp_id" = paste0(chrom, ":", pos)) %>% filter(snp_id %in% prune.set$SNP)

if(nrow(prune.set[!(prune.set$SNP %in% sub.down$snp_id)] == 0))
   {
     message("Mapping complete")
     
}
write.table(sub.down %>% select(snp_id,rsid), file = args[2], quote = FALSE, col.names = FALSE,row.names = FALSE)

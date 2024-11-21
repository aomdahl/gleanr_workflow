library(data.table)
library(dplyr)

full.snps <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/full_variant_qc_metrics.txt.bgz")
full.snps <- full.snps %>% select(chrom, pos, ref, alt, rsid, varid, high_quality, nearest_genes, info, ac_EUR, af_EUR, gnomad_genomes_an_EUR,gnomad_genomes_af_EUR)
sub.snps <- full.snps %>% filter(high_quality == TRUE, nchar(ref) == 1, nchar(alt)==1)
sub.snps$MAF <- sapply(sub.snps$af_EUR, function(x) min(x, 1-x))
sub.snps <- sub.snps %>% filter(MAF > 0.01)
redun <- sub.snps$rsid[duplicated(sub.snps$rsid)]
sub.snps <- sub.snps %>% filter(!(rsid %in% redun)) %>% filter(info > 0.9)
data.table::fwrite(sub.snps, file="/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.gz")

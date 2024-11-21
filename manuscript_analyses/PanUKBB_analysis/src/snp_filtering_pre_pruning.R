#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pacman::p_load(magrittr, dplyr, data.table)
####
#simple script to map SNPs to their rsids from the PanUKBB.
#Designed specifically for the panUKBB
####


#args 1: the input file
#args[1] <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete.union_freq.txt"
#args[2] <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete.1000G_freqs.txt"
snps.lost <- list()

nominees <- fread(args[1], header = FALSE) %>% set_colnames(c("rsid", "freq"))
ref.lookup <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz") 
ref.lookup <- ref.lookup %>% filter(rsid %in% nominees$rsid) %>% 
  select(chrom, pos, ref,alt, rsid, nearest_genes, MAF) %>% 
  mutate("snp_id" =paste0(chrom, ":", pos, ":", ref, ":", alt), "snp_id_inv" = paste0(chrom, ":", pos, ":", alt, ":", ref)) 

message("Number of SNPs lost mapping from nominees to RSIDs (should be 0): ", nrow(nominees) - nrow(ref.lookup))
snps.lost$map_to_rsid <- nrow(nominees) - nrow(ref.lookup)
## A few sanity checks
if(nrow(nominees) > nrow(ref.lookup))
{
  message("Warning: we are unable to map some SNPs in the PanUKBB reference, please be advised.")
}
#Sanity check on the MAF:
if(any(ref.lookup$MAF < 0.01))
{
  message("Warning- SNPs with MAF < 0.01 detected. Possible error in upstream filtering")
}
#Filter down to just hm3 snps
#hm3_ref <- fread("/data/abattle4/aomdahl1/reference_data/hm3_nohla.snpids.bed") %>% set_colnames(c("chr","start","end", "rsid"))
#hm3_ref_alt <- fread("/data/abattle4/aomdahl1/reference_data/hm3_snps_ldsc_ukbb.tsv")
hm3_ref <- fread("/data/abattle4/aomdahl1/reference_data/hapmap_chr_ids.txt")
sub.setted <- (nominees %>% filter(rsid %in% hm3_ref$ID)) #these hm3 ones appear to be enriched for signal
message("Subsetted to HapMap SNPs dropped ", nrow(nominees) -  nrow(sub.setted), " SNPs")
snps.lost$map_to_hm3 <- nrow(nominees) - nrow(sub.setted)
#summary(sub.setted$freq) #worse- loosing some of the shared ones
#summary(nominees$freq) #This is better


nominees.intermediate <- left_join(sub.setted, ref.lookup, by = "rsid") %>% select(snp_id, snp_id_inv, freq, chrom, pos, ref,alt)


#Subset to 1KG SNPs- there its maj:min, so it could be alt:ref or ref:alt
thou.g.ref <- fread("/scratch16/abattle4/ashton/snp_networks/reference_files/thouG_reference_snps.txt", header=FALSE) %>% set_colnames(c("snp_id"))
in.thou.g <- nominees.intermediate %>% filter((snp_id %in% thou.g.ref$snp_id) | (snp_id_inv %in% thou.g.ref$snp_id))
message("Matched to ", nrow(in.thou.g), " SNPs from a query of ", nrow(nominees.intermediate), " in 1KG")
message("Number of SNPs lost mapping from PanUKBB + HM3 to 1KG: ", nrow(nominees.intermediate) - nrow(in.thou.g))
snps.lost$map_to_1KG <- nrow(nominees.intermediate) - nrow(in.thou.g)
#Remove those in HLA region (!)
#coordinates basd on Genome Reference Consortium, https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
message("Filtering out HLA region: coordinates based only hg37, please update if using a different genome build.")
in.thou.g.noHLA <- in.thou.g %>% filter(!(chrom == "6" & pos > 28477796 & pos < 33448354 ))

message("Number of SNPs dropped from the HLA region: ", nrow(in.thou.g) - nrow(in.thou.g.noHLA))
snps.lost$remove_HLA <- nrow(in.thou.g) - nrow(in.thou.g.noHLA)
#report
#Remove ambiguous snps? Not needed here because all of the data was uniformly processed.
message("We are not removing ambiguous SNPs becasue data was processed using a uniform pipeline.")
#If we did remove ambig
ambig.opts <- c("AT", "TA", "CG", "GC")
ambigs = nrow(in.thou.g.noHLA %>% mutate("varvar"=paste0(ref,alt)) %>% filter((varvar %in% ambig.opts)))
#No ambiguous SNPs remain, we have already filtered for those?
if(ambigs == 0)
{
  message("No ambiguous SNPs exist in the data, they have already been filtered out apparently.")
}

#Filtering out multi-allelic SNPs
simple.ids <- paste0(in.thou.g.noHLA$chrom, ":", in.thou.g.noHLA$pos)
dups <-  simple.ids[duplicated(simple.ids)]
in.thou.g.noHLA <- mutate(in.thou.g.noHLA, "temp_snp" = paste0(chrom, ":", pos)) %>% filter(!(temp_snp %in% dups))
message("Number of multi-allelic SNPs dropped: ", length(dups))
snps.lost$multi_allelic <- length(dups)



message("Total number of query SNPs: ", nrow(nominees))
message("Total number of SNPs dropped:",sum(unlist(snps.lost)), " (",round(sum(unlist(snps.lost))/nrow(nominees) * 100, digits = 3), "%)"  )
#This is about 48,000 fewer SNPs lost than in my existing pipeline.


#What about the SNPs 
write.table(in.thou.g.noHLA %>% select(temp_snp, freq) %>% set_colnames(c("SNP","P")) %>% arrange(SNP), 
            file=args[2],
            row.names = FALSE, quote = FALSE)
#"/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete//panUKBB_complete.1000G_freqs.NEW.txt"
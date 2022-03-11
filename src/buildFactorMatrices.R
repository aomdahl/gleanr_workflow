#!/usr/bin/env Rscript
#Manually extract runs from this directory to check for enrichment, see what we get.

writeFactorMatrices <- function(specs1, specs2, studies, run_list, out_prefix)
{
  entry <- c()
  for(i in specs1) { for(j in specs2) { entry <- c(entry, paste0(i, " ", j)) }}
  run_list <- entry
  for(r in run_list)
  {
    index <- which(entry == r)
    curr_run <-  run_stats[[index]]
    f_dat <- data.frame("rownames" = studies, curr_run[[1]])
    #write_tsv(f_dat, paste0(out_prefix, gsub(" ", "_", r), ".factors.txt"))
    write.table(f_dat, file = paste0(out_prefix, gsub(" ", "_", r), ".factors.txt"), sep = "\t", quote = FALSE, row.names=FALSE)
  }
}


#args = commandArgs(trailingOnly=TRUE)
#library(readr)
#out_prefix <- args[2]
#input_data <- args[1]
#specs <- c(scan(text = args[3], what = ""))
#print(specs)
#load(input_data)
#entry <- c()
#studies <- scan("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv", what = character())
#writeFactorMatrices(specs, specs, studies, run_list, out_prefix)

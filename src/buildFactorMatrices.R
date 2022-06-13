#!/usr/bin/env Rscript
#This script has been repursposed as our write_out script. Should be called that, but for convenience staying with this...
#Manually extract runs from this directory to check for enrichment, see what we get.
#s(alphas, lambdas, names, run_stats,output)
writeFactorMatrices <- function(specs1, specs2, studies, snps, run_stats, out_prefix)
{
  entry <- c()
  for(i in specs1) { for(j in specs2) { entry <- c(entry, paste0(i, " ", j)) }}
  run_list <- entry
  for(r in run_list)
  {
    index <- which(entry == r)
    curr_run <-  run_stats[[index]]
    writeFactorMatrix(studies, snps, curr_run, out_prefix)
    
  }
}

writeFactorMatrix <- function(studies, snps, curr_run,identifier, out_prefix)
{
  f_dat <- data.frame("rownames" = studies, curr_run[[1]])
  #write_tsv(f_dat, paste0(out_prefix, gsub(" ", "_", r), ".factors.txt"))
  write.table(f_dat, file = paste0(out_prefix, gsub(" ", "_", identifier), ".factors.txt"), sep = "\t", quote = FALSE, row.names=FALSE)
  #print(paste0(out_prefix, gsub(" ", "_", identifier), ".factors.txt"))
  #how about loadings you ninny?
  if(length(curr_run) > 1)
  {
    if(nrow(curr_run[[2]]) != length(snps))
    {
      message("Unable to complete L matrix, unsure why")
      print(curr_run[[2]])
      #readline()
    } else{
    l_dat <- data.frame("snps" = snps,curr_run[[2]])
    write.csv(x = l_dat, file = paste0(out_prefix, gsub(" ", "_", identifier), ".loadings.txt"), quote = FALSE, row.names=FALSE)
    }

  }else{
    message("No L matrix generated for this run, F too sparse.")
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
#studies <- scan("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv", what = character())
#writeFactorMatrices(specs, specs, studies, run_list, out_prefix)



updateLog <- function(text, verbosity)
{
  if(verbosity == 1)
  {
    log_print(text,console = TRUE)
  } else{
    log_print(text,console = FALSE)
  }
  
}


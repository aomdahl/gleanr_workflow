pacman::p_load(data.table, tidyr, dplyr, magrittr, readr, ggplot2, stringr, optparse)

option_list <- list(
  make_option("--gwas_dir", type = 'character', help = "Directory containing the munged summary stat files", default = ""),
  make_option("--gwas_ext", type = 'character', help = "Extension to the files", default = ".sumstats.gz"),
  make_option("--missing_thresh", type = "numeric", help = "Proportion of NAs per study we will accept", default = 0.1),
  make_option("--output", type = "character", help = "where to write output file and name", default ="missingness.tsv")
)

t=c("--gwas_dir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/saige_benchmark/",
    "--gwas_ext=.sumstats.gz")
#args <- parse_args(OptionParser(option_list=option_list),args=t)
args <- parse_args(OptionParser(option_list=option_list))

inpath <-args$gwas_dir
inext <- args$gwas_ext

files.all <- list.files(path = inpath,pattern = inext )
thresh = args$missing_thresh
missing.count <- c()
any.issues <- FALSE
if(length(files.all) == 0)
{
  message("ERROR: no files for read in")
  print(inpath)
  print(inext)
  quit()
}
for(f in files.all)
{
  test.stats <- fread(paste0(inpath, f))
  #print(paste("Read file", f))
  missingness <- sum(is.na((test.stats$Z)))/nrow(test.stats)
  missing.count <- c(missing.count,missingness )
  if(missingness > thresh)
  {
    message("Warning: missingness exceeds threshold of ", thresh * 100, "%.")
    message("Recommend extracting SNPs directly from summary statistics, not from HM3 subset of SNPs")
    any.issues <- TRUE
  }
}
missingness <- data.frame("file" = files.all, "prop_missing" = missing.count)
if(!any.issues)
{
  message("No files failed missingness test. Congratulations!")
}
write.table(missingness, file = args$output, row.names = FALSE, quote = FALSE, col.names = FALSE)

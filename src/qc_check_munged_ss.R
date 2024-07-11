pacman::p_load(data.table, tidyr, dplyr, magrittr, readr, ggplot2, stringr, optparse)
######
#This script assesses the quality of the outputted summary stats used in downstream analysis, including LDSC.
#This also includes assessing the varinace of the summary stats.
######
option_list <- list(
  make_option("--gwas_dir", type = 'character', help = "Directory containing the munged summary stat files", default = ""),
  make_option("--gwas_ext", type = 'character', help = "Extension to the files", default = ".sumstats.gz"),
  make_option("--missing_thresh", type = "numeric", help = "Proportion of NAs per study we will accept", default = 0.1),
  make_option("--trait_order", type = "character", help = "trait list file, specifying desired order", default = ""),
  make_option("--output", type = "character", help = "where to write output file and name", default ="pheno_report")
)

t=c("--gwas_dir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/",
    "--gwas_ext=.sumstats.gz",
    "--trait_order=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/trait_selections/panUKBB_complete.studies.tsv",
    "--output=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/",
    "--missing_thresh=0.1")
#args <- parse_args(OptionParser(option_list=option_list),args=t)
args <- parse_args(OptionParser(option_list=option_list))

inpath <-args$gwas_dir
inext <- args$gwas_ext
#Another change- just read in what you need to, don't read in everything...
files.all <- list.files(path = inpath,pattern = inext )
thresh = args$missing_thresh
missing.count <- c()
std.dev <- c()
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
  print(paste("Read file", f))
  missingness <- sum(is.na((test.stats$Z)))/nrow(test.stats)
  std.dev <- c(std.dev, sd(test.stats$Z, na.rm = TRUE))
  missing.count <- c(missing.count,missingness )
  if(missingness > thresh)
  {
    message("Warning: missingness exceeds threshold of ", thresh * 100, "%.")
    message("Recommend extracting SNPs directly from summary statistics, not from HM3 subset of SNPs")
    any.issues <- TRUE
  }
}
missingness <- data.frame("file" = files.all, "prop_missing" = missing.count)
sample.var <- data.frame("file" = files.all, "sample_sd" = std.dev)
if(!any.issues)
{
  message("No files failed missingness test. Congratulations!")
}
#match these to the other one...
if(args$trait_order != "")
{
  to = fread(args$trait_order, header = F)
  index.order = c()
  
  #file.names.split = missingness$file %>% gsub(x=., pattern=".sumstats.gz", replacement="") %>% str_split(., pattern = "\\.")
  file.names.split = basename(missingness$file) %>% gsub(x=., pattern=".sumstats.gz", replacement="") %>% 
    str_split(., pattern = "\\.") %>% sapply(., function(x) x[length(x)])
  #The above doesn't work with some ICD codes that have periods in it. More nuanced approach
  for(i in 1:nrow(to))
  {
    pheno = to[i,2]
    match = unlist(lapply(file.names.split, function(x) pheno %in% x))
    if(sum(match) == 1)
    {
      index.order <- c(index.order, which(match))
    }
    else if(sum(match) > 1)
    {
      message("ERROR: files have duplicate names. Output will NOT be in order")
      message("Pipeline may not proceed because 2 studies in the input directory are assigned the same identifying name.")
      print(file.names.split[which(match)])
      #break
      quit()
    }else
    {
      message("No phenotypes match the current file.")
      message("Proceed cautiously, and make sure the names match up.")
      index.order <- c(index.order, i)
    }
  }
  missingness <- missingness[index.order,]
  sample.var <- sample.var[index.order,]
}

write.table(missingness, file = paste0(args$output, "missingness_report.tsv"), row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(sample.var, file = paste0(args$output, "sample_sd_report.tsv"), row.names = FALSE, quote = FALSE, col.names = FALSE)

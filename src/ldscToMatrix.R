pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr,stats, cowplot, argparse, R.utils)
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/postMungeUtils.R")
parser <- ArgumentParser$new()
parser$add_description("Convert LDSC munged summary stats outputs to matrices")
parser$add_argument("--filepathlist", type = "character", help = "Path to 2 column file. First has full file path, 2nd has unique trait name")
parser$add_argument("--feature_list", type = "character", default = "ALL", help = "Path to file specifying just which phenotypes we want. Default is [B,N,SE]")
parser$add_argument("--snp_list", type = "character", help = "Path to list of SNPs to pull out. Default is all those listed in the first file.", default = "")
parser$add_argument("--outdest", type = "character", help = "Path to output destination")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
print(args)

#n.infertility <-getDataFromLDSC("/scratch16/abattle4/ashton/snp_networks/scratch/infertility_analysis/female_infertility_clumped_500kb_r20.05/munged.list.to_trait.txt", 
#                                "N", feature_list="ALL", snp.list = "", fill.nas = TRUE, mean.impute = TRUE)

#z.infertility <-getDataFromLDSC("/scratch16/abattle4/ashton/snp_networks/scratch/infertility_analysis/female_infertility_clumped_500kb_r20.05/munged.list.to_trait.txt", 
                               # "Z", feature_list="ALL", snp.list = NULL, fill.nas = FALSE, mean.impute = FALSE)

if(FALSE)
{
  args <- list()
  args$filepathlist <- "/scratch16/abattle4/ashton/snp_networks/scratch/infertility_analysis/p0.2_FULL/munged.list.to_trait.txt"
  args$feature_list <- "ALL"
  args$snp_list <- ""
}
n.phenos <- getDataFromLDSC(args$filepathlist,"N", feature_list=args$feature_list, snp.list = args$snp_list, fill.nas = TRUE, mean.impute = TRUE)
b.phenos <- getDataFromLDSC(args$filepathlist,"SIGNED_SUMSTAT", feature_list=args$feature_list, snp.list = args$snp_list, fill.nas = FALSE, mean.impute = FALSE)
se.phenos <- getDataFromLDSC(args$filepathlist,"SE", feature_list=args$feature_list, snp.list = args$snp_list, fill.nas = FALSE, mean.impute = FALSE)
message("Writing out....")
write.table(x=n.phenos, file = paste0(args$outdest, "N.tsv"),quote = FALSE, row.names = FALSE)
write.table(x=b.phenos, file = paste0(args$outdest, "B.tsv"),quote = FALSE, row.names = FALSE)
write.table(x=se.phenos, file = paste0(args$outdest, "SE.tsv"),quote = FALSE, row.names = FALSE)

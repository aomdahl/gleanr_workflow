pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr,stats, cowplot, argparse, R.utils)
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/postMungeUtils.R")
parser <- ArgumentParser$new()
parser$add_description("Convert LDSC munged summary stats outputs to matrices")
parser$add_argument("--filepathlist", type = "character", help = "Path to 2 column file. First has full file path, 2nd has unique trait name")
parser$add_argument("--pheno_list", type = "character", default = "ALL", help = "File specifying just which phenotypes we want. Default is ALL")
parser$add_argument("--feature_list", type = "character", default = "ALL", help = "String specifying just which phenotypes we want. Default is ALL, which includes, B,N,SE")
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
features = c(args$feature_list)
if(args$feature_list == "ALL")
{
  c("N", "SE", "SIGNED_SUMSTAT", "Z")
}

for(feature in features)
{
	dat.tab <- getDataFromLDSC(args$filepathlist,feature, feature_list=args$pheno_list, snp.list = args$snp_list, fill.nas = TRUE, mean.impute = TRUE)
	write.table(x=n.phenos, file = paste0(args$outdest, feature, ".tsv"),quote = FALSE, row.names = FALSE)
	
}



pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr,stats, cowplot, optparse, R.utils)
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/postMungeUtils.R")

option_list <- list(
make_option("--filepathlist", type = "character", help = "Path to 2 column file. First has full file path, 2nd has unique trait name"),
make_option("--pheno_list", type = "character", default = "ALL", help = "File specifying just which phenotypes we want. Default is ALL"),
make_option("--feature_list", type = "character", default = "ALL", help = "String specifying just which genetic data we want. Default is ALL, which includes, B,N,SE"),
make_option("--snp_list", type = "character", help = "Path to list of SNPs to pull out. Default is all those listed in the first file.", default = ""),
make_option("--outdest", type = "character", help = "Path to output destination")
)
args <- parse_args(OptionParser(option_list=option_list))

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
features = str_split(args$feature_list, pattern = ",")[[1]] #modify to separate by comma
if(args$feature_list == "ALL")
{
  features = c("N", "SE", "SIGNED_SUMSTAT", "Z")
}

for(feature in features)
{
	dat.tab <- getDataFromLDSC(args$filepathlist,feature, feature_list=args$pheno_list, snp.list = args$snp_list, fill.nas = TRUE, mean.impute = TRUE)
	write.table(x=dat.tab, file = paste0(args$outdest,".", feature, ".tsv"),quote = FALSE, row.names = FALSE)
	
}



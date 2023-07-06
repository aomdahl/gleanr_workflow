#Joins files that were output from LDSC munging..
#I stopped making this because I realized another script already performed this function. Which was it????
#the script you are looking for is ldscToMatrix.R 
library(stringr)
library(data.table)
library(plyr)
library(magrittr)
library(optparse)

option_list <- list(
  make_option(c("--run_list"), type = 'character', help = "Specify a list indicating the File path of each file you wish to join"),
  make_option(c("--col_names"), type = 'character', help = "Specify the corresponding column names for each study, in order. If not given, will default to the file names less the path and extension", default = ""),
  make_option(c("--snps_out"), type = 'character', action = "store_true", help = "Specify which SNPs to extract. If none given, will extract all of them")
)
args <- parse_args(OptionParser(option_list=option_list))

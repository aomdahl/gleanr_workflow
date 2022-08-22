library(stringr)
library(data.table)
dir ="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/"
library(plyr)
library(magrittr)
source(paste0(dir, 'cophenetic_calc.R'))
library(optparse)

option_list <- list(
  make_option(c("--input_path"), type = 'character', help = "Specify the path with all the runs (one directory)"),
  make_option(c("--run_list"), type = 'character', help = "Specify a list indicating the name of each run. Should include the number of K, A, and L"),
  make_option(c("--output"), type = 'character', help = "Specify the name of the ioutput file and full path.", default = "coph.results.txt"),
  make_option(c("--debug"), type = 'logical', action = "store_true", help = "Specify this To run on extisting case.", default = FALSE)
)
args <- parse_args(OptionParser(option_list=option_list))

tail.factors = "_B_SE.1.factors.txt"
tail.loadings = "_B_SE.1.loadings.txt"
if(args$debug)
{
  path = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/infertility/p0.001_FULL/"
  read.files = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/infertility/p0.001_FULL/all.reruns.10iter.txt"

  ofile = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/infertility/p0.001_FULL/coph.summary.10runs.txt"
} else
{
 path = args$input_path; read.files = args$run_list; ofile = args$output
}


#Build the storage structure..
infiles <- fread(read.files,header = FALSE)
#Nested lists- by K, and then by settings
ks <- c()
as <- c()
ls <- c()
for(f in infiles$V1)
{
  ks <- c(ks,str_extract(f, "K\\d+"))
  as <- c(as, str_extract(f, "A\\d+\\.\\d+"))
  ls <- c(ls,str_extract(f, "L\\d+\\.\\d+"))
}
all.matrices <- list()
all.entries <- expand.grid(unique(as), unique(ls))
all.keys <- apply(all.entries, 1, function(x) paste0(x[1], "_", x[2]))
for(k in unique(ks))
{
  all.matrices[[k]] <- list()
  for(e in all.keys)
  {
    all.matrices[[k]][[e]] <- list()
  }
}
for(f in infiles$V1)
{
  k <- str_extract(f, "K\\d+");  as <- str_extract(f, "A\\d+\\.\\d+");  ls <- str_extract(f, "L\\d+\\.\\d+")
  #loadings <- fread(paste0(path, f, tail.loadings))
  factors <- fread(paste0(path, f, tail.factors))
  iter = length( all.matrices[[k]][[paste0(as, "_", ls)]])
  all.matrices[[k]][[paste0(as, "_", ls)]][[(iter + 1)]] <- as.matrix(factors[,-1])
  
}
#Get the counts, make sure we have the thresholds we want.
final.out <- NULL
for(k in unique(ks))
{
  coph.results <- lapply(all.matrices[[k]], function(x) if(length(x) > 1) {copheneticAssessment(x)} else{NA}) #give a list of the 30
  corr.results <- lapply(all.matrices[[k]],  function(x) if(length(x) > 1) {correlationAssessment(x)}else{NA}) #give a list of the 30
  corr.results.per.k <- lapply(all.matrices[[k]],  function(x) if(length(x) > 1) {scaledCorrelationAssessment(x)}else{NA}) #give a list of the 30
  nruns <- sapply(all.matrices[[k]], function(x) if(length(x) > 0) {length(x)}else{NA})
  stopifnot(names(coph.results) == names(corr.results))
  
  final.out <- rbind(final.out, data.frame("Coph" = unlist(coph.results), "Corr_frob" = unlist(corr.results),"Corr_per_factor"=unlist(corr.results.per.k), "Settings" = names(coph.results), "K" = k, "nruns" = unlist(nruns)))
  
}
write.table(final.out,file = ofile, quote = FALSE, row.names = FALSE)
#recommend the best one
library(dplyr)
library(tidyr)
top <- (data.frame(final.out) %>%drop_na() %>% filter(Coph > 0.9) %>% arrange(Corr_frob))[1,]
message("Top recommended settings")
print(top)

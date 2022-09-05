library(stringr)
library(data.table)
dir ="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/"
library(plyr)
library(magrittr)
source(paste0(dir, 'cophenetic_calc.R'))
source(paste0(dir, 'sparsity_scaler.R'))
library(optparse)

#HElper function- storage structure for analysis:
organizeData <- function(as, ls)
{
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
  return(all.matrices)
}
#lengthData(all.matrices.f)
lengthData <- function(ds)
{
  tl <- 0
  for(k in 1:length(ds))
  {
    for(param in 1:length(ds))
    {
      tl <- tl + length(ds[[k]][[param]])
    }
  }
  tl
}
#readIndata(infiles$V1, as, ls, tail.factors)
readIndata <- function(infiles, as, ls, tail.str="_B_SE.*.factors.txt")
{
  read.files <- c()
    #read in the 
  data.structure <- organizeData(as, ls)
  total.read.in <- 0
  for(f in infiles)
  {
    k <- str_extract(f, "K\\d+");  as <- str_extract(f, "A\\d+\\.\\d+");  ls <- str_extract(f, "L\\d+\\.\\d+")
    matching.files <- Sys.glob(paste0(path, f,tail.str), dirmark = FALSE)
    for(mf in matching.files)
    {
      if(!(mf %in% read.files))
      {
        factors <- fread(mf)
        iter = length(data.structure[[k]][[paste0(as, "_", ls)]])
        data.structure[[k]][[paste0(as, "_", ls)]][[(iter + 1)]] <- as.matrix(factors[,-1])
        read.files <- c(read.files, mf)
        total.read.in <- total.read.in + 1
      }else
      {
        message("Skipping file, already read in")
        mf
      }

    }
  }
  message("Read in a total of: ", total.read.in, " files.")
  return(data.structure)

}




option_list <- list(
  make_option(c("--input_path"), type = 'character', help = "Specify the path with all the runs (one directory)"),
  make_option(c("--run_list"), type = 'character', help = "Specify a list indicating the name of each run. Should include the number of K, A, and L"),
  make_option(c("--output"), type = 'character', help = "Specify the name of the ioutput file and full path.", default = "coph.results.txt"),
  make_option(c("--debug"), type = 'logical', action = "store_true", help = "Specify this To run on extisting case.", default = FALSE)
)
args <- parse_args(OptionParser(option_list=option_list))

tail.factors = "_B_SE.*.factors.txt"
tail.loadings = "_B_SE.*.loadings.txt"
if(args$debug)
{
  path = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/infertility/p0.001_FULL/subsampling_1/"
  read.files = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/infertility/p0.001_FULL/subsampling_1/run_list.txt"

  ofile = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/infertility/p0.001_FULL/coph.results.avgk17.UPDATED.txt"
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
#Get the full parameter space.
for(f in infiles$V1)
{
  ks <- c(ks,str_extract(f, "K\\d+"))
  as <- c(as, str_extract(f, "A\\d+\\.\\d+"))
  ls <- c(ls,str_extract(f, "L\\d+\\.\\d+"))
}
message("Reading in F matrices...")
all.matrices.f <- readIndata(infiles$V1, as, ls, tail.factors)
#message("Read in a total of, " lengthData(all.matrices.f))
message("Reading in L matrices...")
all.matrices.l <- readIndata(infiles$V1, as, ls, tail.loadings)

#Get the counts, make sure we have the thresholds we want.
message("Assessing performance across all...")
final.out <- NULL
for(k in unique(ks))
{
  coph.results <- lapply(all.matrices.f[[k]], function(x) if(length(x) > 1) {copheneticAssessment(x)} else{NA}) #give a list of the 30
  corr.results <- lapply(all.matrices.f[[k]],  function(x) if(length(x) > 0) {correlationAssessment(x)}else{NA}) #give a list of the 30
  corr.results.per.k <- lapply(all.matrices.f[[k]],  function(x) if(length(x) > 0) {scaledCorrelationAssessment(x)}else{NA}) #give a list of the 30
  nruns <- sapply(all.matrices.f[[k]], function(x) if(length(x) > 0) {length(x)}else{NA})
  #non-empty entries

  non.empty.count <- sapply(1:length(all.matrices.f[[k]]), function(i) if(length(all.matrices.f[[k]][[i]]) > 0) {nonEmptyAvg(all.matrices.f[[k]][[i]],all.matrices.l[[k]][[i]])} else{{NA}})
  #matrix sparsity.
  f.sparsity <- lapply(all.matrices.f[[k]], function(x) if(length(x) > 0) {matrixSparsityAvg(x)} else{NA}) #give a list of the 30
  l.sparsity <- lapply(all.matrices.l[[k]], function(x) if(length(x) > 0) {matrixSparsityAvg(x)} else{NA}) #give a list of the 30
  stopifnot(names(coph.results) == names(corr.results))
  final.out <- rbind(final.out, 
    data.frame("Coph" = unlist(coph.results), "Corr_frob" = unlist(corr.results),
      "Corr_per_factor_pair"=unlist(corr.results.per.k), "Settings" = names(coph.results), "Init_K" = k, "Actual_K" = non.empty.count, 
      "F_sparsity" = unlist(f.sparsity), "L_sparsity" = unlist(l.sparsity), "nruns" = unlist(nruns)))
  
}
write.table(final.out,file = ofile, quote = FALSE, row.names = FALSE)
#recommend the best one
library(dplyr)
library(tidyr)
top <- (data.frame(final.out) %>%drop_na() %>% filter(Coph > 0.9) %>% arrange(Corr_frob))[1,]
message("Top recommended settings")
print(top)

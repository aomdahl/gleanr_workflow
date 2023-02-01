#!/usr/bin/Rscript
library(stringr)
library(data.table)
dir ="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/"
#library(plyr)
library(dplyr)
library(magrittr)
source(paste0(dir, 'cophenetic_calc.R'))
source(paste0(dir, 'sparsity_scaler.R'))
source(paste0(dir, 'read_in_tools.R'))
source(paste0(dir, 'pve.R'))
source(paste0(dir, 'l1_parameter_selection.R'))
library(optparse)

readInDataForEval <- function(args)
{
  bargs <- args
  bargs$covar_matrix <- ""
  bargs$genomic_correction <- ""
  bargs$weighting_scheme <- "B"
  bargs$gwas_effects <- args$betas
  bargs$trait_names <- ""
  bargs$IRNT <- FALSE
  readInData(bargs)
}



#This is copy/pased from the BIC code, with some important modifications noted.

quickSort <- function(tab, col = 1)
{
  tab[order(tab[,..col], decreasing = TRUE),]
}
#This will get the pve across the list of matrices.
avgPVE <- function(full.dat, m)
{
  mean(sapply(m, function(i) sum(PercentVarEx(full.dat, as.matrix(i)))))
}

#we don't want to be redundant, so omit the ones that are repeated...
#alphas: the optimal parameters to consisder
#all.a: the full list of parameters from which they are drawn
#n.points: how many new points to pick from.
ProposeSparsityParamsFromGrid <- function(alphas,all.a, n.points)
{
  #N-points will include the original ones...
  if(length(alphas) == 1)
  {
    new.list <- singlePointGridExpansion(all.a[order(all.a)], alphas, n.points)
  }else
  {
    range.a <- range(alphas)
    new.list <- (seq(range.a[1],range.a[2],length = n.points))
  }
  if(any(is.na(new.list)))
  {
    message("detected possible ERROR in proposed parameters")
    print(new.list)
    message("Dropping NAs")
    new.list <- new.list[!is.na(new.list)]
  }
  if(any(new.list < 0))
  {
    message('Detected negatives in proposed parameters; dropping these')
    new.list <- new.list[new.list > 0]
  }
  if(length(new.list) == 0)
  {
    message("Possible error, no terms detected")
    message("Retruning original input.")
    return(alphas)
  }
return(new.list)

}


ProposeNewSparsityParamsFromScore <- function(bic.list,sparsity.params, n.points = 7, no.score = FALSE)
{
  #we want to keep it griddy. So a bit further downstream compare this to the BIC, I think this actually the one that we want, but okay.
  if(no.score)
  {
    #then bic.list is the optimal one; generate fake scores
    fake.scores <- rep(100,length(sparsity.params))
    fake.scores[which(sparsity.params == bic.list)] <- -1
    bic.list <- fake.scores
  }
  optimal.sparsity.param <- sparsity.params[which.min(bic.list)]
  message("optimal sparsity param is ", optimal.sparsity.param)
  ordered.list <- sparsity.params[order(sparsity.params)] #order the list
  #New paradigm: always look above and below,
  #If its the smallest paramter tested
  singlePointGridExpansion(sparsity.params[order(sparsity.params)], optimal.sparsity.param, n.points)
}

#HElper function- storage structure for analysis:
organizeData <- function(as, ls, ks)
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
#all.matrices.f <- readIndata(infiles$V1, as, ls,path, tail.factors)
#readIndata(infiles$V1, as, ls, tail.factors)

readIndata <- function(infiles, as, ls,ks, path, tail.str="_B_SE.*.factors.txt")
{
  read.files <- c()
    #read in the 
  data.structure <- organizeData(as, ls, ks)
  total.read.in <- 0
  for(f in infiles)
  {
    k <- str_extract(f, "K\\d+");  as <- str_extract(f, "A\\d+\\.*\\d+");  ls <- str_extract(f, "L\\d+\\.*\\d+")
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
      }
      #else
      #{
      #  #message("Skipping file, already read in")
      #  mf
      #}

    }
  }
  message("Read in a total of: ", total.read.in, " files.")
  return(data.structure)

}

#we are running the script.
if (sys.nframe() == 0){
 
#TODO: modify so look at cophenetic correlation across all those with a given K
option_list <- list(
  make_option(c("--input_path"), type = 'character', help = "Specify the path with all the runs (one directory)"),
  make_option(c("--run_list"), type = 'character', help = "Specify a list indicating the name of each run. Should include the number of K, A, and L"),
  make_option(c("--output"), type = 'character', help = "Specify the name of the ioutput file and full path.", default = "coph.results.txt"),
  make_option(c("--debug"), type = 'logical', action = "store_true", help = "Specify this To run on extisting case.", default = FALSE),
  make_option(c("--betas"), type = 'character', help = "Specify original betas. Needed to get PVE", default = NULL),
  make_option(c("--scale_n"), type = 'character', help = "Specify Sample size to scale by.", default = NULL),
  make_option(c("--lambda_gc"), type = 'character', help = "Specify lambda gc. Optional to get PVE", default = NULL)
)

args <- parse_args(OptionParser(option_list=option_list))

tail.factors = "_B_SE.*.factors.txt"
tail.loadings = "_B_SE.*.loadings.txt"
#also a test of the subsampling method, huh.
if(args$debug)
{
  args <- list()
  #="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/renv_f/results/ukbb_GWAS_h2-0.1_rg-0.9/model_selection_nov2022_renv/"
  #path="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/renv_f/results/ukbb_GWAS_h2-0.1_rg-0.9/model_selection_nov2022_log_grid/"
  path = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/udler_original/model_grid_search_k6/"
  #path="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/renv_f/results/ukbb_GWAS_h2-0.1_rg-0.9/model_selection_nov2022_renv/2nd_iter/"
  #read.files = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_GWAS_h2-0.1_rg-0.9/factorization/full_list_night_oct18.txt"
  #read.files = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//results/ukbb_GWAS_h2-0.1_rg-0.9/model_selection_nov2022/runs.so.far.txt"
  #read.files="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/renv_f/results/ukbb_GWAS_h2-0.1_rg-0.9/model_selection_nov2022_log_grid/all.runs.dec5.txt"
  read.files="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/udler_original/model_grid_search_k6/full.run.list.txt"
  
  #read.files = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/renv_f/results/ukbb_GWAS_h2-0.1_rg-0.9/model_selection_nov2022_renv/2nd_iter/all.param.runs.txt"
  #ofile = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/infertility/p0.001_FULL/coph.results.oct_18.txt"
  ofile = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/udler_original/model_grid_search_k6/run_performance.txt"
  #args$betas <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/B.tsv"
  args$betas <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/beta_signed_matrix.tsv"
  args$scale_n <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/processed_data/sample_counts_matrix.tsv"
  #args$lambda_gc <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_GWAS_h2-0.1_rg-0.9/Lambda_gc.tsv"
} else
{
 path = args$input_path; read.files = args$run_list; ofile = args$output
}

all <- readInDataForEval(args)
X <- all$X * all$W

#Build the storage structure..
infiles <- fread(read.files,header = FALSE)
#print(infiles)
#Nested lists- by K, and then by settings
ks <- c()
as <- c()
ls <- c()
#Get the full parameter space.
for(f in infiles$V1)
{
  ks <- c(ks,str_extract(f, "K\\d+"))
  as <- c(as, str_extract(f, "_A[\\.0-9e-]+_"))%>% gsub(., pattern = "_", replacement = "")
  ls <- c(ls,str_extract(f, "_L[\\.0-9e-]+$")) %>% gsub(., pattern = "_", replacement = "")
}
message("Reading in F matrices...")
files <- paste0(path, infiles$V1)
all.matrices.f <- readIndata(infiles$V1, as, ls,ks, path, tail.str = tail.factors)
#message("Read in a total of, " lengthData(all.matrices.f))
message("Reading in L matrices...")
all.matrices.l <- readIndata(infiles$V1, as, ls,ks, path, tail.str = tail.loadings)
a.n <- unique(as.numeric(gsub(as, pattern = "A", replacement = "")))
l.n <- unique(as.numeric(gsub(ls, pattern = "L", replacement = "")))
#Get the counts, make sure we have the thresholds we want.
message("Assessing performance across all...")
final.out <- NULL
true.k.groups <- list()
group.tracking <- list()
message("Getting key statistics about the runs...")
for(k in unique(ks))
{
  k.n <- as.numeric(gsub(k, pattern = "K", replacement = ""))
  if(k.n == 0)
  {
    message('Init K unknown, going with first matrix present')
    k.n <- ncol(all.matrices.f[[k]][[1]][[1]])
  }
  coph.results <- lapply(all.matrices.f[[k]], function(x) if(length(x) > 1) {copheneticAssessment(x)} else{NA}) #need more than 1 matrix to do this....
  corr.results <- lapply(all.matrices.f[[k]],  function(x) if(length(x) > 0) {correlationAssessment(x)}else{NA}) #give a list of the 30
  corr.results.per.k2 <- lapply(all.matrices.f[[k]],  function(x) if(length(x) > 0) {scaledCorrelationAssessment(x, type = "squared")}else{NA})
  #The median of the averages
  corr.avg <- lapply(all.matrices.f[[k]],  function(x) if(length(x) > 0) {scaledCorrelationAssessment(x, type = "average_r2")}else{NA})
  corr.results.per.k <- lapply(all.matrices.f[[k]],  function(x) if(length(x) > 0) {scaledCorrelationAssessment(x, type = "std")}else{NA})#give a list of the 30
  nruns <- sapply(all.matrices.f[[k]], function(x) if(length(x) > 0) {length(x)}else{NA})
  #non-empty entries

  non.empty.count <- sapply(1:length(all.matrices.f[[k]]), function(i) if(length(all.matrices.f[[k]][[i]]) > 0) {nonEmptyAvg(all.matrices.f[[k]][[i]],all.matrices.l[[k]][[i]])} else{{NA}})
  #true.numFactors <- sapply(1:length(all.matrices.f[[k]]), function(i) if(length(all.matrices.f[[k]][[i]]) > 0) {nonEmptyAvg(all.matrices.f[[k]][[i]],all.matrices.l[[k]][[i]])} else{{NA}})
  #matrix sparsity.
  f.sparsity <- lapply(all.matrices.f[[k]], function(x) if(length(x) > 0) {matrixSparsityAvg(x,k.n )} else{NA}) #give a list of the 30
  l.sparsity <- lapply(all.matrices.l[[k]], function(x) if(length(x) > 0) {matrixSparsityAvg(x,k.n)} else{NA}) #give a list of the 30
  #PVE:
    pve <- lapply(all.matrices.f[[k]], function(m) if(length(m) > 0 & !is.null(args$betas)) {avgPVE(X, m)} else{NA})
  
  stopifnot(names(coph.results) == names(corr.results))
  final.out <- rbind(final.out, 
    data.frame("Coph" = unlist(coph.results), "Corr_frob" = unlist(corr.results),
      "Corr_per_factor_pair"=unlist(corr.results.per.k),"Corr_per_factor_pair_SQ"=unlist(corr.results.per.k2), "Med_Average_R2" = unlist(corr.avg), 
      "Settings" = names(coph.results), "Init_K" = k, "Group_nonempty_average" = non.empty.count, 
      "F_sparsity" = unlist(f.sparsity), "L_sparsity" = unlist(l.sparsity), "nruns" = unlist(nruns), "Avg_pve"=unlist(pve)))
  #need the cohenetic correlation across all from the same group
  for(mats in all.matrices.f[[k]])
  {
    for(m in mats)
    {
      k.true <- ncol(m)
      if(!(k.true %in% group.tracking))
      {
        group.tracking <- c(group.tracking, k.true)
        true.k.groups[[k.true]] <- list("1"=m)
      }
      else
      {
        true.k.groups[[k.true]][[length( true.k.groups[[k.true]]) + 1]] <- m
      }
    }
  }
}

#Look at the cophenetic correlation per group of true
coph.by.k <- lapply(true.k.groups, function(x) if(length(x) > 0) {copheneticAssessment(x)} else{NA})
n.by.k <- lapply(true.k.groups, function(x) length(x))
m <- data.frame("Actual_K" = 1:length(coph.by.k), "K_group_coph" = unlist(coph.by.k), "N" = unlist(n.by.k))
#The true stability of each matrix at a given K is importna
write.table(m %>% filter(!is.na(K_group_coph)),file = paste0(ofile, "k_coph.txt"), quote = FALSE, row.names = FALSE)
#12/20: ADD IN UPDATED VERSION OF WHAT I THINK YUAN WAS DOING... ....
final.out <- final.out %>% tidyr::separate(Settings, into = c("Alpha", "Lambda"), sep= "_",remove = FALSE) #%>% left_join(., m, by = "Actual_K")
by.group.med <- final.out %>% group_by(Group_nonempty_average) %>% summarize("Coph_med" = median(Coph))
final.out <- left_join(final.out, by.group.med, by = "Group_nonempty_average")
write.table(final.out,file = paste0(ofile, "all.txt"), quote = FALSE, row.names = FALSE)

#######3recommend the best one
#As of Dec 2022, conditions are:
# K.
library(tidyr)
#global K needs to be > 0.9
# K' > 0.9
# F sparsity needs to be > 0.1 and < 0.8
# PVE needs to be > 20%?
# Of these options, choose the one that has minimal off-diagonal R2
#acceptable.k <- (m %>% filter(K_group_coph > 0.9))$Actual_K
#message("K settings that pass threshold:")
#print(paste(acceptable.k))
#tops <- (data.frame(final.out) %>% tidyr::drop_na() %>% filter(Coph > 0.9) %>% filter(F_sparsity > 0.3, F_sparsity < 0.9) %>% 
#          arrange(Med_Average_R2) %>% filter(Group_nonempty_average %in% acceptable.k, Med_Average_R2 < 0.01, Avg_pve > 0.5, L_sparsity > 0.05))
#adopting the more lenient interpretation of what I believe Yuan was doing...
#another change- don't put a hard R2 filter in there, just pick the minimal ones (?)
tops <- (data.frame(final.out) %>% tidyr::drop_na() %>% filter(Coph > 0.9) %>% filter(F_sparsity > 0.3, F_sparsity < 0.9) %>% 
          arrange(Med_Average_R2) %>% filter(!is.na(Coph_med), Coph_med > 0.9, Avg_pve > 0.5, L_sparsity > 0.05))
if(nrow(tops) == 0)
{
  acceptable.k <- (m %>% filter(K_group_coph > 0.8))$Actual_K
  tops <- (data.frame(final.out) %>% tidyr::drop_na() %>% filter(Coph > 0.9) %>% filter(F_sparsity > 0.3, F_sparsity < 0.9) %>% 
             arrange(Med_Average_R2) %>% filter(Group_nonempty_average %in% acceptable.k, Avg_pve > 0.5, L_sparsity > 0.05))
}

message("Top recommended settings")
print(tops)
#in the case of multiple settings that come to the top, these give a decent range.
message("If proceeding to next step, next settings are given by.....")
best.a <- unique(as.numeric(gsub(tops$Alpha, pattern = "A", replacement = "")))
best.l <- unique(as.numeric(gsub(tops$Lambda, pattern = "L", replacement = "")))
if(nrow(tops) >= 1)
{
  new.a <- ProposeSparsityParamsFromGrid(best.a,a.n, 4)
  new.l <- ProposeSparsityParamsFromGrid(best.l,l.n,4)
  print("Alpha:")
  print(paste0((new.a),collapse = ","))
  print("Lambda:")
  print(paste0((new.l),collapse = ","))
}else #(nrow(tops) == 0)
{
  message("paramters too strict, select manually")
}
#visualize approach
visualize = FALSE
if(visualize)
{
  library(ggplot2)
  old.grid <- expand.grid(a.n, l.n) %>% mutate("s" = "initial grid")
  new.grid <- expand.grid(new.a, new.l)%>% mutate("s" = "new grid")
  df <- data.frame(rbind(new.grid, old.grid)) %>% arrange(Var1, Var2) %>% 
    mutate("s" = ifelse(Var2 %in% best.l & Var1 %in% best.a, "optimal", s)) %>% unique()
  #First grid:
  ggplot(df %>% filter(s %in% c("initial grid", "optimal")), aes(x = Var1, y= Var2, color = s)) + geom_point(size = 3) + 
    theme_minimal() + xlab(bquote(alpha)) + ylab(bquote(lambda)) + 
    labs(color = "Grid search Iteration") + ggtitle("1st Grid search")
  
  ggplot(df, aes(x = Var1, y= Var2, color = s)) + geom_point(size = 3) + 
    theme_minimal() + xlab(bquote(alpha)) + ylab(bquote(lambda)) + 
    labs(color = "Grid search Iteration") + ggtitle("Grid search visualization")
  
  
  plot(tops$F_sparsity, ends$L_sparsity, xlim = c(0,1), ylim = c(0,1), 
       pch = 19, xlab = "Sparsity F", ylab = "Sparsity L")
  plot(tops$Avg_pve, ends$Med_Average_R2, xlim = c(0,1), ylim = c(0,0.1), 
       pch = 19, xlab = "PVE", ylab = "Factor correlation")
  
}
}


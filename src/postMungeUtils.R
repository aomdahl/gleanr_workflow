
#Take munged LDSC summary stats and get a list of some column out from them
#l is the list of files and their names
#col.out is the variable you want out (either Z,P, or N usually)
#feature_list is the name of the traits you want, shold match with column 2 of entries in l
#snp.list is list of RSIDs to keep.
#returns entliust
library(data.table)
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)

#rename from fullnanes
getPhenoName <- function(stin){
  unlist(str_split(string = stin, pattern = "\\.") %>% lapply(., function(x) x[[(length(x)-2)]]))
}


#Picks the element largest in list 1; if it corresponds with multiple entries, pick the ones that is smallest in list 2
maxNminSE <- function(x1, x2)
{
  mn <- max(x1)
  if(length(which(x1 == mn)) > 1)
  {
    opt <- x2[which(x1 == mn)]
    r <- min(opt)
    return(which(x2 == r)[which(x2 == r) %in% which(x1 == mn)])
  }
  return(which.max(x1))
}

#Goal of this is to get data across all of the different runs for each SNP
#Reports the distribution statistics of the paramter of interest and the missingness rate

snpDistsFromLDSC <- function(l, col.out, feature_list="ALL", snp.list = "", save="")
{
  
  #l <- "/scratch16/abattle4/ashton/snp_networks/scratch/infertility_analysis/p0.2_FULL/munged.list.to_trait.txt"
  l <- "/scratch16/abattle4/ashton/snp_networks/scratch/infertility_analysis/hm3_extracts/munged.list.to.trait.txt"
  d <- listReadInHelper(l, feature_list, snp.list)
  nums <- list()
  lin <- d$file.list
  snp.list <- unlist(d$snps)
  feature_list = d$traits
  snp.stats <- matrix(NA, nrow = length(snp.list), ncol =length(feature_list) )
  for(i in 1:nrow(lin))
  {
    curr <- fread(lin$V1[i]) %>% select(SNP, !!sym(col.out))
    print(i)
    if(all(curr$SNP == snp.list))
    {
      snp.stats[,i] <- unlist(curr[,2])
    }
    else{

        message("SNPs not lined up...")
    }
    
  }
  #now, get the summary stats on each one....
  r <- data.frame("SNPs" =snp.list, "NA_count"=apply(snp.stats, 1, function(x) sum(is.na(x))), 
             "Max" =apply(snp.stats, 1, function(x) max(x,na.rm = TRUE)), "Min" = apply(snp.stats, 1,function(x) min(x,na.rm = TRUE)), 
             "Avg" = rowMeans(snp.stats, na.rm = TRUE),  "Var" = apply(snp.stats, 1,function(x) var(x,na.rm = TRUE)))
  #if 
  #save = "./scratch16-abattle4/ashton/snp_networks/scratch/infertility_analysis/hm3_extracts/z_score_collective.stats.tsv"
  if(save != ""){write.table(r,out , quote = FALSE,row.names = FALSE)}
#plot(r$Max, r$NA_count)
  return(r)

}

listReadInHelper <- function(l, feature_list = "ALL", snp.list = "")
{
  lin <- fread(l,header = FALSE)
  if(feature_list != "ALL")
  {
    lin <- lin %>% filter(V2 %in% feature_list) %>% mutate("ordered" = factor(V2, levels = feature_list)) %>% arrange(ordered)
  }else
  {
    feature_list = unlist(lin[,2])
  }
  if(snp.list == "")
  {
    message("No SNPs given- going by first entry in the list")
    snp.list <- fread(lin$V1[1])$SNP
  }
  
  snp.list <- data.frame("SNP" = snp.list)
  return(list("snps" = snp.list, "file.list" = lin, "traits" = feature_list))
}

getDataFromLDSC <- function(l, col.out, feature_list="ALL", snp.list = "", fill.nas = FALSE, mean.impute = TRUE)
{
 #l <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/pairwise.traits.labels.txt"
  #l <- "/scratch16/abattle4/ashton/snp_networks/scratch/infertility_analysis/female_infertility_clumped_500kb_r20.05/munged.list.to_trait.txt"

  d <- listReadInHelper(l, feature_list, snp.list)
  nums <- list()
  lin <- d$file.list
  snp.list <- d$snps
  feature_list = d$traits
  
  nums <- lapply(1:nrow(lin), function(x) fread(lin$V1[x]) %>% left_join(snp.list, ., by = "SNP") %>% select(!!sym(col.out)))
  tot <- as.matrix(do.call("cbind", nums))
  #Mean impute the sample size for SNPs with no specified N.....
  if(fill.nas)
  {
    cms <- rep(0, ncol(tot))
    if(mean.impute)
    {
      cms <- colMeans(tot, na.rm = T)
    }
    
    #Replace the NAs with the appropriate value
    for(i in 1:ncol(tot))
    {
      tot[,i][is.na(tot[,i])] <- cms[i]
    }
    stopifnot(!any(is.na(tot)))
  }
  tot <- cbind(snp.list, tot)
  colnames(tot) <- c("SNP", feature_list)

  data.frame(tot)
  #write.table(tot, "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/variant_lists/hm3.pruned.N.txt", quote = FALSE, row.names = FALSE)
}


#pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr,stats, cowplot)
snp.list <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/variant_lists/hm3.pruned.N.txt")$SNP
##Tool for processing the results from LDSC for weights....
#This requires having run pairwise LDSC already
#ldsc.path- where the ldsc output files are from my pipeline
#file.extension is what teh path ending for those ldsc things look like
#n.snps- how many snps are in each file? that we want to replicate?


lambdaGCWeights <- function(ldsc.path, file.list, snp.list, feature_list,  file.extension = "*ldsc_report.csv")
{
  n.snps <- length(snp.list)
  #file.list <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/pairwise.traits.labels.txt"
  lin <- fread(file.list,header = FALSE) %>%
    filter(V2 %in% feature_list) %>% mutate("ordered" = factor(V2, levels = feature_list)) %>% arrange(ordered) %>% 
    mutate("study" = basename(V1))
  #ldsc.path = "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/ldsr_results/"
  ldsc.main.files <- list.files(ldsc.path, pattern = file.extension)
  joined.main <- do.call("rbind", lapply(ldsc.main.files, 
                                         function(x) fread(paste0(ldsc.path, "/", x)) %>% mutate("Source" = x)))
  joined.main$study <- gsub( joined.main$study, pattern = ".bgz", replacement = ".gz")
  joined.main.max.snp <- joined.main %>% group_by(study) %>% slice(which.max(snp_overlap)) %>% 
    mutate("upper" = intercept + intercept_se * 1.96,
                                 "lower" = intercept - intercept_se * 1.96) %>% 
    mutate("contains_1" = ifelse(lower < 1 & upper > 1, TRUE, FALSE)) %>% 
    mutate("intercept_use" = ifelse(contains_1 | is.na(intercept) | intercept < 1 , 1, intercept)) %>%
    select(study, intercept_use) %>% left_join(., lin, by = "study") %>% arrange(ordered) %>% ungroup() %>%
    select(ordered, intercept_use) %>% drop_na()
  
  correction.matrix <- do.call("rbind", lapply(1:n.snps, function(x) joined.main.max.snp$intercept_use))
  colnames(correction.matrix) <- joined.main.max.snp$ordered
  
  write.table(data.frame("SNP" = snp.list, correction.matrix), 
              "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/variant_lists/hm3.pruned.GC.txt", quote = FALSE, row.names = FALSE)

}


ldscTableToMatrix <- function(joined_main, col_pick, diag_scores = 1, null_scores = 0)
{
  overlap.mat <- joined_main %>% select(p1, p2, !!sym(col_pick)) %>% mutate("id" = paste0(p1, ":", p2)) %>% distinct_at(vars(id), .keep_all = TRUE) %>% 
    select(-id) %>% pivot_wider(values_from = !!sym(col_pick), names_from =p1) %>% arrange(p2) 
  
  #which ons are missing?
  #colnames(overlap.mat)[!(colnames(overlap.mat) %in% overlap.mat$p2)]
  #overlap.mat$p2[!(overlap.mat$p2 %in% colnames(overlap.mat))]
  #Squareify it.
  new_mat <- list()
  names_list <- c()
  for(i in 2:ncol(overlap.mat))
  {
    curr_frick <- colnames(overlap.mat)[i]
    names_list <- c(names_list, curr_frick)
    if(curr_frick %in% overlap.mat$p2)
    {
      new_mat[[i-1]] <- overlap.mat[which(overlap.mat$p2 == curr_frick),]
    }else{
      print("missing this...")
      print(curr_frick)
      new_mat[[i-1]] <- rep(NA, ncol(overlap.mat))
    }
  }
  #Make this into a matrix, then add in the missing ones as a new row and nuew column
  new.overlap.mat <- data.frame(do.call("rbind", new_mat))
  colnames(new.overlap.mat) <- colnames(overlap.mat)
  missing.entries <- which(!(overlap.mat$p2 %in% colnames(overlap.mat)))
  for (m in missing.entries)
  {
    new.overlap.mat <- rbind(new.overlap.mat, overlap.mat[m,])
    str.name <- unlist(overlap.mat[m,1]$p2)
    new.overlap.mat <- cbind(new.overlap.mat, c(t(overlap.mat[m,-1]), 1))
    last.entry <- ncol(new.overlap.mat)
    colnames(new.overlap.mat)[last.entry] <- str.name
  }
  stopifnot(all(colnames(new.overlap.mat)[-1] == new.overlap.mat$p2))
  ldscintmatrix <- as.matrix(new.overlap.mat[,-1])
  diag(ldscintmatrix) <- diag_scores

  #make sure symmetric- garabage code
  m <- ldscintmatrix
  for(i in 1:nrow(m))
  {
    for(j in 1:ncol(m))
    {
      if( j >= i) {
        if(is.na(m[i,j]) |is.na(m[j,i]) |  (round(m[i,j], digits = 3) != round(m[j,i], digits = 3)))
        {
          #message("non-symmetric entry. Replacing if NA, leaving otherwise.")
          if(is.na(m[i,j]) & (!is.na(m[j,i])))
          {
            m[i,j] <- m[j,i]
          } else if(is.na(m[j,i]) & !is.na(m[i,j]))
          {
            m[j,i] <- m[i,j]
          }
        }else
        {
          #message("No replacement possible for cell.")
        }
      }
    }
  }
  m[is.na(m)] <- null_scores
  return(m)
}


quickRead <- function(look_path, pattern)
{
  #look_path <- "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/ldsr_results"
  ldsc.data <- list.files(look_path,pattern = pattern)
  #expanding this for saftey....
  #inter <- lapply(ldsc.data, function(x) fread(paste0(look_path, x)) %>% mutate("Source" = x))
  inter <- list()
  i=1
  for(x in ldsc.data)
  {
    file = fread(paste0(look_path, x))
    if(nrow(file) > 0)
    {
      inter[[i]] <- file %>% mutate("Source" = x)
      i=i+1
    }
    else
    {
      message("An empty table passed in for, ", look_path)
    }
  }
  do.call("rbind", inter)
}
#look_path = args$input_path
#which.col = "h2_obs"
#Does this line up rows and columns correctly?
#ldscGCOVIntoMatrix(look_path = args$input_path,"rg")
ldscGCOVIntoMatrix <- function(look_path = "/scratch16/abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/ldsr_results/", which.col, 
                               filter_se = FALSE, filter_fdr = 1)
{
  joined.main <- quickRead(look_path, pattern = "*rg_report.tsv") %>% mutate("Source2" = gsub(Source, pattern = ".rg_report.tsv", replacement = ""))
  joined.count.info <- quickRead(look_path, pattern = "*.ldsc_report.csv")%>% mutate("Source2" = gsub(Source, pattern = ".ldsc_report.csv", replacement = ""))
  joined.main$study <- basename(joined.main$p2)
  joined.main$lookup <- paste0(joined.main$study, ":", joined.main$Source2)
  joined.count.info$lookup <- paste0(joined.count.info$study, ":", joined.count.info$Source2)
  all.dat <- left_join(joined.main, joined.count.info, by = "lookup") %>% 
    select(p1, p2, rg, se, z, p, h2_obs, h2_obs_se, h2_int, h2_int_se, gcov_int, gcov_int_se,Source.x, snp_overlap) %>% 
    rename("Source" = Source.x)
  joined.main <- all.dat
  
  if(which.col == "rg_se")
  {
    #ldscTableToMatrix(joined.main,"se", null_scores = 1)
    ns = 0
    which.col = "se"
    
  }
  
  
  if(filter_se)
  {
    message("Nullifying entries that are within the 95% CI")

    
    if(which.col == "gcov_int") #is the intercept 0 or not?
    {
      #note- because we have a large sample size, and our df is n-1, we can use the normal approximation for the t
      #
      if(filter_fdr != 1)
      {
	      message("estimating a global FDR for correction...")
      z.stat <- joined.main$gcov_int/joined.main$gcov_int_se
      pvals.sub <- 2*pnorm(abs(z.stat[!is.na(z.stat)]), lower.tail = F)
      fdr.sub <- p.adjust(pvals.sub, method = "fdr")
      fdr <- rep(1, length(z.stat))
      fdr[!is.na(z.stat)] <- fdr.sub #put in the FDR.
      joined.main$gcov_int[fdr > filter_fdr] <- 0 #keep only those passing my FDR estimate.
      } else(filter_se)
      {
        #Option to do it with confidence interval, or with a p-value of some kind of correction.
        upper <- joined.main$gcov_int + 1.96*joined.main$gcov_int_se;
        lower <- joined.main$gcov_int - 1.96*joined.main$gcov_int_se
        nullify_cases <- (lower < 0 & upper > 0) | is.na(joined.main$gcov_int)
        joined.main$gcov_int[nullify_cases] <- 0
      }

      
    }
    else
    {
      message("no filtering to be done on specified column")
    }
  }
  
  if(which.col == "h2_obs" | which.col == "h2_int")
  {
    #Choosing the ones with the largest sample sizes and smallest SE.
    #tmp <- joined.main %>% filter(h2_obs != "NA" & !is.na(h2_obs)) %>% group_by(p2) %>% slice(maxNminSE(snp_overlap,h2_obs_se)) %>% 
    #  select(-gcov_int_se, -Source, -gcov_int, -snp_overlap) %>% ungroup()
    tmp <- joined.main %>% filter(h2_obs != "NA" & !is.na(h2_obs)) %>% group_by(p2) %>% arrange(h2_obs_se) %>% slice(which.max(snp_overlap)) %>% 
      select(-gcov_int_se, -Source, -gcov_int, -snp_overlap) %>% ungroup()
    #To get N

    
    renames <- unlist(str_split(string = tmp$p2, pattern = "\\.") %>% lapply(., function(x) x[[(length(x)-2)]]))
    tmp$Trait = renames
    return(tmp %>% ungroup() %>% select(-p1))
  }
  #If null scores...
  ns = 0
  if(which.col == "p")
  {
    ns = 1
  }
  ldscTableToMatrix(joined.main,which.col, null_scores = ns)
}

#which.col doesn't matter.
#ldscStdIntoMatrix(look_path = args$input_path,"" , filter_se = args$filter_se)
ldscStdIntoMatrix <- function(look_path, which.col, filter_se = FALSE)
{
  #look_path="/scratch16/abattle4/ashton/snp_networks/scratch/infertility_analysis/hm3_ldsc/tabular/"
	print(look_path)
  joined.main <- quickRead(look_path, pattern = "*.ldsc_report.csv" ) %>% filter(!is.na(intercept))
  #sanity checck on the numbers- some longer because ran after the fact..
  #For each study, choose teh one with the largest sample size on the RUNS.
  #question- does the # of overlapping snps correspond to the 
  #Note: if the run is using the full set of SNPS of each time, then the estimates should be identical from run to run.
  #quick test: are the stats identical for a trait? NO
  # the standard error DROPS for  SummaryStatistics_E2_WOMEN_META_191025.estradiol_2.sumstats.gz as N increases
  #stick with guns- largest N snp overlap is good.
  ci <- 1 
  #get the one with the most samples
  #get the intercept if the standard error doesn't capture 1
  # fix the names
  max.n <- joined.main %>% group_by(study) %>% slice(which.max(snp_overlap)) %>% 
    mutate("ci_bands" = ifelse((((intercept + ci*intercept_se) > 1) & ((intercept - ci*intercept_se) < 1)) | (intercept < 1), 1, intercept)) %>%
    ungroup()
  max.n$phenotype <- getPhenoName(max.n$study)

  #make sure all our SNP overlaps are good
  #If the number of SNPs is below 200,000, this is usually bad, and ldsc will print a warning. This is what they say on the wiki
  if(any(max.n$snp_overlap < 200000))
  {
    message("Beware- studies have fewer than 200,000 SNPS contributing to their estimates. This is recommended against by LDSC")
  }
  if(filter_se)
  {
    #If +- 1 standard error includes 1, juust make it 1. Otherwise, let it be.
    ret.tab <- select(max.n, phenotype, ci_bands) 
  }else
  {
    ret.tab <- select(max.n, phenotype,intercept )
  }
  
  return(ret.tab %>% set_colnames(c("phenotype", "intercept_keep")))
}






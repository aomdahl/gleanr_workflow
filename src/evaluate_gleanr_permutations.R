#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization//src/factorEvaluationHelper.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/getTopGenes_helper.R")
###############


#' Wrapper for the logistic regression enrichment tests
#'
#' @param u_vector simple vector of effect sizes from U
#' @param u_gene_ids corresponding (same order) genes for each entry of U
#' @param gene_set the list of genes in the set we want to test for
#'
#' @return a single row of a datframe
#' @export
#'
#' @examples
logit_enrichment_test <- function(u_vector,u_gene_ids, gene_set) {
  if (length(gene_set) != length(unique(gene_set))) {
    warning("Gene set entries aren't unique. Could be a problem with input.")
  }
  max.logit <- testEnrichmentLogit(gene_set,u_gene_ids, abs(u_vector), collapse = "max")
  mean.logit <- testEnrichmentLogit(gene_set,u_gene_ids, abs(u_vector), collapse = "mean")
  
  data.frame("maxlogit_p"= summary(max.logit$test)$coef[8],"maxlogit_b"=max.logit$logOR,"maxlogit_SE"=max.logit$logOR_SE,
             "avglogit_p"=summary(mean.logit$test)$coef[8],"avglogit_b"=mean.logit$logOR,"avglogit_SE"=mean.logit$logOR_SE)
}
#  

# enrichment_test(gene_set,set_id, matrices_list[[i]], matrix_v2g, num_col, bg_genes="all")
#target_gene_set=gene_set
#target_gene_set_id=set_id
#matrix=matrices_list[[i]]
#factor_index=num_col
#bg_genes="all"
#' Performs all the target 
#'
#' @param target_gene_set - biological set (e.g. GO terms) you want to test for
#' @param target_gene_set_id - an identifier label to name that set
#' @param matrix - the U matrix you are testing on
#' @param matrix_v2g - ordered list of RSIDs and genes (names(rsid and gene_id) that match U)
#' @param factor_index - the column of U to evaluate
#' @param bg_genes - background set of genes
#'
#' @return
#' @export
#'
#' @examples
enrichment_test <- function(target_gene_set,target_gene_set_id, matrix, matrix_v2g, factor_index,bg_genes="all",...)
{
  if (length(target_gene_set) != length(unique(target_gene_set))) {
    warning("Gene set entries aren't unique. Could be a problem with input.")
  }
  
  #Perform continuous test with logistic regression
  logit_res <- logit_enrichment_test(u_vector=matrix[,factor_index], u_gene_ids=matrix_v2g$gene_id, gene_set=target_gene_set)
  
  #Perform fisher's exact test:
  ## First get the top snps, genes, and background snps, genes
  top.snps <- prioritize_snps(as.matrix(matrix[,factor_index]),snp_ids = matrix_v2g$rsid,method="top_fe", plot_curve = FALSE)
  #top_factor_genes <- map_snps_to_genes(top.snps[[1]], matrix_v2g) #i think maybe this is wrong
  #TODO- try to use the one above, has some error catching built in.
  top_factor_genes = (matrix_v2g %>% dplyr::filter(rsid %in% top.snps[[1]]))$gene_id
  if(bg_genes == "all") {
    bg_genes = unique(matrix_v2g$gene_id)
  }else if(bg_genes == "bottom"){
    bg_genes = unique(matrix_v2g$gene_id[!(matrix_v2g$gene_id %in% top_factor_genes)])
  }else{
    message("Using a pre-specified bg list you supplied.")
  }
  
  fisher_res  <- testEnrichment(unique(top_factor_genes), unique(bg_genes), target_gene_set,conf.level=0.9,...)
  fisher_res_pseudo  <- testEnrichment(unique(top_factor_genes), unique(bg_genes), target_gene_set,conf.level=0.9,pseudocount=TRUE,...)
  #Return in a single row of a data frame.
  cbind(data.frame("factor" = factor_index,"gene_set"=target_gene_set_id, "fisher_p"=fisher_res$test$p.value,"fisher_OR"=fisher_res$test$estimate,
                   "fisher_CI_upper" = fisher_res$test$conf.int[2], "fisher_CI_lower" = fisher_res$test$conf.int[1], 
                   "fisher_logOR"=fisher_res$logOR, "fisher_logOR_SE"=fisher_res$logOR_SE,
                   "pseudofisher_p"=fisher_res_pseudo$test$p.value,"pseudofisher_OR"=fisher_res_pseudo$test$estimate,
                   "pseudofisher_logOR"=fisher_res_pseudo$logOR, "pseudofisher_logOR_SE"=fisher_res_pseudo$logOR_SE,
                   "pseudofisher_CI_upper" = fisher_res_pseudo$test$conf.int[2], "pseudofisher_CI_lower" = fisher_res_pseudo$test$conf.int[1]),
        logit_res)
}



##############
# Define command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to RData file containing list of matrices", metavar = "FILE"),
  make_option(c("-g", "--gene_sets"), type = "character", help = "Path to file with gene sets", metavar = "FILE"),
  make_option(c("-n", "--num_columns"), type = "character", help = "Comma-separated list of numbers of columns to test", metavar = "STRING"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory", metavar = "DIR"),
  make_option(c("-r", "--ref_analysis"), type = "character", help = "Path to reference analysis", metavar = "FILE"),
  make_option(c("-t", "--ref_analysis_only"), type = "logical", help = "SPecify this if you want to end after running oon the real data. ",default=FALSE, action="store_true")
)

t =c("--input=/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed23_/permute_100.RData",
  "--gene_sets=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/manuscript_analyses/PanUKBB_analysis/analysis/platelet_diseases_gene_lists.txt",
  "--num_columns=4,11,23","--output=/scratch16/abattle4/ashton/snp_networks/scratch/manuscript_reviews/permute_testing/assess_test1",
  "--ref_analysis=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData")

bd="/scratch16/abattle4/ashton/snp_networks/"
t2 = c(paste0("--input=", bd, "/scratch/manuscript_reviews/permute_testing/panUKBB_41K_final/seed27_byROW/permute_200.RData"),
       paste0("--gene_sets=", bd, "/custom_l1_factorization/manuscript_analyses/PanUKBB_analysis/analysis/platelet_diseases_gene_lists.txt"),
       "--num_columns=4,11,23",
       paste0("--output=", bd, "/scratch/manuscript_reviews/permute_testing/assess_test_seed27_200_byROW/"),
       paste0("--ref_analysis=", bd, "//custom_l1_factorization/results/panUKBB_complete_41K_final/panUKBB_complete_41K_final_final_dat.RData"))

opt <- parse_args(OptionParser(option_list = option_list))
#opt <- parse_args(OptionParser(option_list=option_list), args = t2)
# Validate required arguments
if (is.null(opt$input) || is.null(opt$gene_sets) || is.null(opt$num_columns) || is.null(opt$output) || is.null(opt$ref_analysis)) {
  stop("All arguments --input, --gene_sets, --num_columns, --ref_analysis, and --output are required.", call. = FALSE)
}

# Load input data
load(opt$input)
if (!exists("shuff_results_U")) stop("Expected 'shuff_results_U' object in input file.")
matrices_list <- shuff_results_U

#load reference data
load(opt$ref_analysis)
if (!exists("ret")) stop("Expected 'ret' object in input file.")

#Load gene sets
gene_sets_file <- fread(opt$gene_sets,header = FALSE)
gene_sets <- lapply(1:nrow(gene_sets_file), function(i) unlist(strsplit(gene_sets_file$V2[i], ",")))
names(gene_sets) <- gene_sets_file$V1

if(opt$num_columns == "all")
{
  num_columns <- as.numeric(1:ncol(matrices_list[[1]]))
}else
{
  num_columns <- as.numeric(unlist(strsplit(opt$num_columns, ",")))
}


#make the gene set 
matrix_v2g <- load_v2G() %>% mutate("rsid" = factor(rsid, levels = ret$snp.ids)) %>% dplyr::filter(rsid %in% ret$snp.ids) %>% arrange(rsid) %>% dplyr::select(rsid, gene_id)
library(dplyr)


#Optional: write out the results of the "true" test set
message("First performing statistical tests on the reference data...")
true.results <- list()
for (num_col in num_columns) {
  for(j in seq_along(gene_sets)) {
    set_id=names(gene_sets)[j]
    gene_set = gene_sets[[j]]
    true.results[[paste0("Cols", num_col, "_GeneSet",j)]] <- 
      enrichment_test(gene_set,set_id,ret$U, matrix_v2g, num_col, bg_genes="all") 
    #Some stuff for debugging
    if(num_col == 4 & j==5)
    {
      f4_set.new <- unique(gene_set)
      u.dat <- ret$U
      or_test_hmt <-  true.results[[paste0("Cols", num_col, "_GeneSet",j)]]
      bg_genes="all"
      save(f4_set.new,u.dat, matrix_v2g,num_col,or_test_hmt,bg_genes="all", file="/scratch16/abattle4/ashton/snp_networks/scratch/review_responses/enrichment_test_hmtp_NEW_debug.RData")
      
    }
  }
}
output_file_true <- file.path(opt$output, "enrichment_results.tsv")
write.table(do.call(rbind, true.results), file = output_file_true, quote = FALSE, row.names = FALSE, sep = "\t")
message("Tested the reference data.")
if(opt$ref_analysis_only)
{
  message("Ending now.")
  quit()
}

# Run permuted enrichment tests
results <- list()
for (i in seq_along(matrices_list)) {
  for (num_col in num_columns) {
    for(j in seq_along(gene_sets)) {
      if(is.null( matrices_list[[i]])) #this occurs b/c the 2nd 200 get stored in indices 101-200, so the first 100 list entries are empty.
      {
        if(i==1)
        {
          message("Skipping empty results...")
        }
       
      }else
      {
        set_id=names(gene_sets)[j]
        gene_set = gene_sets[[j]]
        #message(paste0("Matrix", i, "_Cols", num_col, "_GeneSet",j))
        results[[paste0("Matrix", i, "_Cols", num_col, "_GeneSet",j)]] <- 
          enrichment_test(gene_set,set_id, matrices_list[[i]], matrix_v2g, num_col, bg_genes="all") 
      }

    }
  }
}

# Save results
output_file <- file.path(opt$output, "permuted_enrichment_results.tsv")
write.table(do.call(rbind, results), file = output_file, quote = FALSE, row.names = FALSE, sep = "\t")

message("Enrichment analysis complete. Results saved to ", output_file)

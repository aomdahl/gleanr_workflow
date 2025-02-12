# Load necessary libraries
pacman::p_load(tidyverse, optparse, data.table, magrittr, EnsDb.Hsapiens.v79)
#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/getTopGenes_helper.R")
source(normalizePath("src/getTopGenes_helper.R"))


#opt$snp_gene_map,opt$snp_id_map, opt$snp_scores, opt$method, opt$output
#snp_gene_map_file <- opt$snp_gene_map; snp_id_map <- opt$snp_id_map; snp_scores_file <- opt$snp_scores; method <- opt$method; odir <- opt$output
main <- function(snp_gene_map_file,snp_id_map, snp_scores_file, method, odir) {
  #FIRST STEPS= MAP SNPS TO GENES
  # Read input files
  message("Loading in reference data (includes openTargets distances)")
  message(snp_gene_map_file)
  snp_gene_map <- data.frame(loadSNPtoGeneMap(snp_gene_map_file))
  
  #Upgraded version of joined singular.. a little more straightforward
  print(colnames(snp_gene_map))
  snp_gene_map.singular <- snp_gene_map %>% group_by(hg38) %>% slice_min(n=1,order_by = d) %>% ungroup() %>%
    dplyr::filter(nchar(ref_allele) == 1, nchar(alt_allele)==1) %>% distinct(across(c(hg38,gene_id,overall_score)), .keep_all=TRUE)
  
  #any remaining duplicates or missing genes- either NA for distance or distances exactly the same (weird)
  multi.snps <- snp_gene_map.singular %>% group_by(hg38) %>% summarize("count"=n()) %>% dplyr::filter(count>1)
  missing.genes <- snp_gene_map.singular %>% dplyr::filter(is.na(gene_id))
  snps.to.fix <- unique(c(multi.snps$hg38, missing.genes$hg38))
  tss_coords = "/data/abattle4/aomdahl1/reference_data/gencode/gencode.v38.basic.annotation.gene_transcripts.TSS.bed" #Just assumet his is right for now, circle back later.
  top.for.reduns <- tiebreakers_by_distance(dplyr::filter(snp_gene_map.singular, hg38 %in%snps.to.fix),tss_coords)
  snp_gene_map.singular <- rbind(snp_gene_map.singular %>% dplyr::filter(!(hg38 %in% snps.to.fix)), top.for.reduns)
  stopifnot(length(unique(snp_gene_map.singular$hg38)) == length(snp_gene_map.singular$hg38))
  
  #Still have 40 of these. Match by distance if possible, otherwise just pick the first one.
  #singular.mapping <- loadAltSNPtoGeneMap(bed_path="/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/snp_to_gene/gencode.v38.basic.annotation.closest_tss.bed",snp_gene_map)
  snp_id_map.df <- fread(snp_id_map) %>% mutate("hg38"=  gsub(x=paste0(V1,":",V3), pattern="chr", replacement = ""))
 
  
  #joint with the RSID
  joined.dat <- left_join(snp_gene_map, snp_id_map.df) #Used in write out only, keeps repeats
  joined.singular <- left_join(snp_gene_map.singular, snp_id_map.df,by="hg38") %>%
    set_colnames(c("hg38", "ref","alt","overall_score","d", "gene_id","snp_chr","snp_start","snp_end", "rsid")) %>%
    dplyr::filter(!is.na(rsid))
  snp_scores <- fread(snp_scores_file)
  #snp_scores <- snp_scores %>% dplyr::filter(!SNPs %in% missing.snps)
  snp_ids <- unlist(snp_scores[,1])
  missing.snps <- snp_ids[which(!(snp_ids %in% joined.singular$rsid))] #lost via liftover
  hg37.missed.snps <- recoverUnmappedSNPs(missing.snps) #assuming the extract here is correct, take some time to verify

 
  snp_scores <- as.matrix(snp_scores[,-1])
  #Have a nice, clean, combined source of all the SNPs that is unique:
   clean.joined <- rbind(joined.singular %>% dplyr::select(hg38, rsid, ref, alt, overall_score, gene_id),
                         hg37.missed.snps %>% dplyr::select(-distance))
   #Save a copy for the supplement and for other reference:
   write_ref="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/gene_snp_mappings.csv"
   if(!file.exists(write_ref))
   {
     write.out.joined <- rbind(joined.singular %>% dplyr::select(hg38, rsid, ref, alt, overall_score, gene_id,d) %>% dplyr::rename("distance"=d),
                               hg37.missed.snps) %>% dplyr::select(rsid, hg38, ref, alt, gene_id, overall_score, distance) %>%
       dplyr::rename("OpenTargets_score"=overall_score)
     write.table(write.out.joined,file="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/gene_snp_mappings.csv", sep=",",row.names = FALSE,
                 quote = FALSE)
   }

  stopifnot(length(unique(clean.joined$rsid)) == nrow(clean.joined))
  stopifnot(nrow(clean.joined) == nrow(snp_scores))
  
    # Prioritize SNPs
    message("Prioritizing snps now....")
    prioritized_snps <- prioritize_snps(snp_scores, snp_ids, method)
      #bg_snps <- lapply(prioritized_snps, function(x) snp_ids[!(snp_ids %in% x)])
    bg_snps <- lapply(prioritized_snps, function(x) snp_ids) #should be all snps, the full possible set.
    
    # Map prioritized SNPs to genes
    # function(prioritized_snps, snp_gene_map, snp_id_map)
    #prioritized_genes <- map_snps_to_genes(prioritized_snps, singular.mapping,snp_id_map.df) #erm....
    prioritized_genes <- map_snps_to_genes(prioritized_snps, clean.joined) #erm....
    bg_genes <- map_snps_to_genes(bg_snps, clean.joined)
    message("Writing output....")
    #write it out for each factor
    write_tabular_reports(snp_scores,snp_ids, prioritized_genes, prioritized_snps,joined.dat,clean.joined, bg_genes, odir)


  
}

# Define command-line arguments
option_list <- list(
  make_option(c("-g", "--snp_gene_map"), type = "character",
              help = "CSV file mapping SNPs to genes", metavar = "FILE", 
              default="/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/snp_to_gene/41K_openTargets.withDist.csv"),
  make_option(c("-s", "--snp_scores"), type = "character", default = NULL, 
              help = "Factor U file containing SNP scores matrix"),
  make_option(c("-m", "--method"), type = "character", default = "factor_elbow", 
              help = "SNP prioritization method (global_elbow, factor_elbow, top_percent, top_fe)", metavar = "METHOD"),
  make_option(c("--snp_id_map"), type = "character",  
              help = "File mapping hg38 SNP ids to RSIDs. Should be output from a previous run of liftover.", 
              default="/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/snp_to_gene/local_liftover.hg38.sorted.bed"),
  make_option(c("--snp_gene_map_alt",type="character",
                help="BED file with alternative SNP gene mapping if desired",
                default="/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/snp_to_gene/bedtools_nearest_coding_match.bed")),
  make_option(c("--output"), type = "character",  
              help = "Path for output files")
  
)

test = c("--snp_scores=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/latent.loadings.txt",
         "--method=top_fe", "--output=./here")
# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
#opt <- parse_args(opt_parser, args = test)
opt <- parse_args(opt_parser)
# Check for required arguments
if (is.null(opt$snp_gene_map) || is.null(opt$snp_scores) || is.null(opt$method)) {
  print_help(opt_parser)
  stop("Please provide all required arguments", call. = FALSE)
}

# Run main script
main(opt$snp_gene_map,opt$snp_id_map, opt$snp_scores, opt$method, opt$output)






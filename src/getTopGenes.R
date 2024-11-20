# Load necessary libraries
pacman::p_load(tidyverse, optparse, data.table, magrittr, EnsDb.Hsapiens.v79)
##FUNCTIONS FOR PRIORITINZG SNPS
# Function to check matches
checkMatch <- function(bic.dat) {
  sapply(bic.dat$rec.dat$Vs, function(x) {
    if(ncol(x) == ncol(bic.dat$optimal.v)) {
      return(all(x == bic.dat$optimal.v))
    } else {
      return(FALSE)
    }
  })
}

# Function to check distances
checkDist <- function(bic.dat) {
  sapply(bic.dat$rec.dat$Vs, function(x) {
    if(ncol(x) == ncol(bic.dat$optimal.v)) {
      return(sum((x - bic.dat$optimal.v)^2))
    } else {
      return(-1)
    }
  })
}

# Function to plot weight distribution, based on the Udler approach
plotWeightDist <- function(input.mat, title, nsample=10000, log10=FALSE) {
  vals <- as.vector(as.matrix(input.mat))
  if(length(vals) > nsample) {
    message("Downsampling vals for visualization")
    vals.samp <- sample(vals, size = nsample)
  } else {
    vals.samp <- vals
  }
  png("test.png")
  dev.control(displaylist="enable") 
  if(log10) {
    plot(rank(-vals.samp), log10(vals.samp+1), pch = 19, col = "skyblue", xlab = "Weight number", ylab = "Clustering weight", main = title)
  } else {
    plot(rank(-vals.samp), vals.samp, pch = 19, col = "skyblue", xlab = "Weight number", ylab = "Clustering weight", main = title)
  }
  p <- recordPlot()
  dev.off()
  return(p) 
}

# Function to calculate deltas, based on the udler approach
calcDeltas <- function(input.mat, sort=TRUE) {
  vals <- as.vector(as.matrix(input.mat))
  d <- length(vals)
  ordered <- vals[order(-vals)]
  if(sort) {
    sort(sapply(2:d, function(x) ordered[x -1] - ordered[x]), decreasing = TRUE)
  } else {
    (sapply(2:d, function(x) ordered[x -1] - ordered[x]))
  }
}

# Function to plot deltas, based on the udler approach
plotDeltas <- function(input.mat, plot.title) {
  deltas <- calcDeltas(input.mat)
  d <- length(deltas)
  png("test.png")
  dev.control(displaylist="enable")  
  plot(1:(d-1), deltas, pch = 19, col = "skyblue", xlab = "Rank" , ylab = "Change in clustering weight", main = plot.title)
  p <- recordPlot()
  dev.off()
  return(list("plot"=p, "deltas"=deltas))
}

# Function to find cutoff score, based on the udler approach
findCutoffScore <- function(input.mat, perc.thresh = 0.05) {
  vals <- as.vector(as.matrix(input.mat))
  d <- length(vals)
  ordered <- sort(vals, decreasing = TRUE)
  deltas <- sapply(2:d, function(x) ordered[x -1] - ordered[x])
  delta.cutoff <- quantile(deltas, probs = c(1-0.05))
  min_thresh <- max(which(deltas > delta.cutoff))
  message("last significant delta: ", delta.cutoff)
  ordered[min_thresh]
}

# Function to get vector elbow
getVectorElbow <- function(v, plot_curve =TRUE, prop.keep=1) {
  sorted.valsdf <- data.frame("x"=rank(-v), "y"=v) %>% arrange(x)
  elbow.i <- ceiling(pathviewr::find_curve_elbow(sorted.valsdf, plot_curve = plot_curve) * prop.keep)
  sorted.valsdf[elbow.i, 2 ]
}
getTopPercent <- function(v, perc = 0.01)
{
  quantile(v,probs=c(1-perc))
}

# Define a function to prioritize SNPs based on the selected method
#prioritize_snps(snp_scores, snp_ids, method)
prioritize_snps <- function(snp_scores,snp_ids, method) {
  # Placeholder for SNP prioritization logic
  if (method == "global_elbow") {
    # Implement global_elbow prioritization logic
    global.elbow <- getVectorElbow(abs(c(snp_scores)))
    top.snp.indices <- apply(snp_scores, 2, function(x) which(abs(x) > global.elbow))
    message("not TESTED")
  } else if (method == "factor_elbow") {
    # Implement factor_elbow prioritization logic
    all.elbows <- apply(snp_scores, 2, function(x) getVectorElbow(abs(x)))
    top.snp.indices <- lapply(1:length(all.elbows), function(i) which(abs(snp_scores[,i]) > all.elbows[i]))
  } else if (method == "top_percent") {
    # Implement percent prioritization logic
   all.top.cutoff <- apply(snp_scores, 2, function(x) getTopPercent(abs(x)))
    top.snp.indices <- lapply(1:length(all.top.cutoff ), function(i) which(abs(snp_scores[,i]) > all.top.cutoff [i]))
  } else if(method == "top_fe"){
    #Get upper third of elbow genes
    all.elbows.half <- apply(snp_scores, 2, function(x) getVectorElbow(abs(x), prop.keep = 1/3))
    
    top.snp.indices <- lapply(1:length(all.elbows.half), function(i) which(abs(snp_scores[,i]) > all.elbows.half[i]))
  }
    else {
    stop("Invalid prioritization method")
  }
  
 prioritized_snps <- lapply(top.snp.indices, function(x) snp_ids[x])

  return(prioritized_snps)
}

################################
############
# FUNCTIONS FOR MAPPING SNPS TO GENES
############
################################
# Define a function to map SNPs to genes
#prioritized_genes <- map_snps_to_genes(prioritized_snps, singular.mapping,snp_id_map.df) #erm....
#clean_joined
#(bg_snps, clean.joined)

map_snps_to_genes <- function(prioritized_snps, snp_gene_map)#snp_gene_map, snp_id_map) {
{
  #Report overall # missing SNPs
  # Placeholder for SNP to gene mapping logic
  #gotten.snps <- lapply(prioritized_snps, function(x) (snp_id_map %>% dplyr::filter(V4 %in% x))$hg38) #map the rsids to hg38. BOO.
  #Loosing some SNPs here. This would be due to lost liftover. Not many though
  #lost.snps <- sapply(1:length(prioritized_snps), function(i) abs(length(gotten.snps[[i]]) - length(prioritized_snps[[i]])))
  #if(any((lost.snps) > 0))
  #{
  #  message("As many as ", max((lost.snps)), " lost when lifting-over onto hg38 for gene sets")
  #}
  #gene_mappings <- lapply(gotten.snps, function(x) (snp_gene_map %>% dplyr::filter( hg38 %in% x) %>% 
  #                                                   group_by(hg38) %>% slice_head(n=1) %>% ungroup() %>%
  #                                                    dplyr::filter(!is.na(gene_id)))$gene_id)
  #NEW VERSION
  gene_mappings <- lapply(prioritized_snps, function(x) (snp_gene_map %>% dplyr::filter(rsid %in% x))$gene_id)
  stopifnot(sum(sapply(gene_mappings, function(x) sum(is.na(x)))) == 0)

  hist(sapply(gene_mappings, length), xlab="# Genes", main="Histograms of # genes per factor",breaks = 20) #do this as a bar chart

  #mapping incorrect $
  lost.snps <- sapply(1:length(gene_mappings), function(i) abs(length(gene_mappings[[i]]) - length(prioritized_snps[[i]])))
  prop.lost <- lost.snps / sapply(prioritized_snps, length)
  message("Max percent lost snps in gene sets: ", max(prop.lost))
  return(gene_mappings)
}


loadSNPtoGeneMap <- function(path)
{
  snp_gene_map <- fread(path)
  
  # 1. Convert from ensembl.gene to gene.symbol
  
  if(grepl(pattern="ENSG", x=snp_gene_map$gene_id[1]))
  {
    message("converting Ensemble to gene symbols")
    conversion.genes <- (unique(snp_gene_map$gene_id))
    geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= conversion.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID")) %>% 
      dplyr::rename("gene_id"=GENEID)
    snp_gene_map<- left_join(snp_gene_map, geneIDs1, by= "gene_id") %>% dplyr::select(-gene_id) %>% dplyr::rename("gene_id"=SYMBOL)
  }
  #Now pick the highest scoring one per SNP!
  snp_gene_map %>% group_by(hg38) %>% slice_max(order_by=overall_score)
}

#shoot.... some SNPs don't get mapped at all, waht about them?
#loadAltSNPtoGeneMap(bed_path="/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/snp_to_gene/gencode.v38.basic.annotation.closest_tss.bed",snp_gene_map)
#DEPRECATED, DON'T USE THIS
loadAltSNPtoGeneMap <-function(bed_path, curr_map)
{
  
  tss_coords = "/data/abattle4/aomdahl1/reference_data/gencode.v38.basic.annotation.gene_transcripts.TSS.bed" #Just assumet his is right for now, circle back later.
  #Loaded data mapping each SNP to its closest gene
  bed.map <- fread(bed_path,header = FALSE) %>% set_colnames(c("snp_chr","snp_start","snp_end","rsid", "gene_chr","gene_start","gene_end", "gene", "distance")) %>% 
    mutate("hg38" = paste0(gsub(x=snp_chr,pattern = "chr",replacement = ""),":",snp_end)) %>% distinct()
  #There are 89 SNPs where the distances are the same. remember that.
  multis <- bed.map %>% group_by(rsid) %>% summarize("count"=n()) %>% dplyr::filter(count >1)
  #First, map SNPs that have no current gene map to the nearest one
  no.match <- curr_map %>% dplyr::filter(is.na(gene_id)) %>% group_by(hg38) %>% slice_head(n=1) %>% ungroup() %>% left_join(., bed.map, by = "hg38") %>% 
    dplyr::select(hg38,ref_allele,alt_allele,overall_score, gene,rsid, distance) %>% distinct()
  no.match$overall_score <- NA
  
  #Now, for the duplicates, where opentargets maps one snp to multiple genes
  multi.snps <- (curr_map %>% dplyr::filter(!is.na(gene_id)) %>% group_by(hg38) %>% summarize("count"=n()) %>% dplyr::filter(count >1))$hg38
  tss.coords <- fread(tss_coords) %>% set_colnames(c("chr","start","end","gene")) %>% 
    mutate("hg38" = paste0(gsub(x=chr,pattern = "chr",replacement = ""),":",end))
  
  #Extract just the SNPs mapping to multiple
  #I don't like doing it this way.... too many options to mess up...
  multi.snp.dat <-   curr_map %>% dplyr::filter(hg38 %in% multi.snps) %>% separate(hg38,into = c("chr","pos"),sep=":",remove = FALSE)
  for(snp in multi.snps)
  {
    
    cd <- multi.snp.dat %>% dplyr::filter(hg38 %in% snp)
    stopifnot(length(unique(cd$overall_score))==1)
    chr_ <- paste0("chr", unique(cd$chr))
    pos_ = as.numeric(unique(cd$pos))
    top_hits <-(tss.coords %>% dplyr::filter(chr==chr_) %>% mutate("distance"=abs(end-pos_)) %>% arrange(distance))
    matched.genes <- dplyr::filter(top_hits, gene %in% cd$gene_id)
    if(nrow(matched.genes) == 0)
    {
      message("No gene matches, just going with the top 1....")
      ret.hit <- cd[1,]
    }else
    {
      ret.hit <- cd %>% dplyr::filter(gene_id == matched.genes$gene[1])
    }
    #SNPs with multiple maps
    
    #For each group, pick the
    remaining.issues <- test %>% dplyr::filter(matches_candidate == "none") %>% dplyr::select(hg38,overall_score, gene_id, rsid, snp_chr,snp_end) %>% 
      group_by(hg38) %>% group_modify(~mapGeneMinDist(.))
    
    #Okay, at this point, all SNPs have been mapped,where possible. eat it.
    #Now, need to get everything to the right format
    #hg38, ref_allele, alt_allele, overall_score, gene_id
    #Create vector of solo maps
    
    fixed.mappings <- rbind(no.match %>% rename(gene="gene_id") %>% dplyr::select(hg38, ref_allele, alt_allele, overall_score, gene_id),
                            solved %>%mutate("ref_allele" = NA, "alt_allele"=NA) %>% dplyr::select(hg38, ref_allele, alt_allele, overall_score, gene_id),
                            remaining.issues   %>%mutate("ref_allele" = NA, "alt_allele"=NA) %>% dplyr::select(hg38, ref_allele, alt_allele, overall_score, gene_id)) %>%
      distinct() %>% group_by(hg38) %>% slice_max(overall_score)
    fine.mappings <- curr_map %>% dplyr::filter(!(hg38 %in% fixed.mappings$hg38))
    return(rbind(fine.mappings, fixed.mappings))
  }  

  
  
  test <- curr_map %>% dplyr::filter(hg38 %in% multi.snps) %>% 
    left_join(.,bed.map %>% dplyr::select(hg38,rsid,gene_chr,gene_start,gene_end,gene,snp_chr, snp_end,distance),by="hg38")  %>% group_by(hg38) %>% slice_max(overall_score) %>%
    mutate("matches_candidate" = ifelse(any(gene_id == gene), "match", "none")) %>% ungroup()
  
  #Those cases where we have a gene that matches open targets, choose the mapping where the distances is the minimum.
  solved <- test %>% dplyr::filter(matches_candidate == "match") %>% dplyr::select(-ref_allele, -alt_allele) %>% dplyr::filter(gene_id == gene) %>% distinct()%>%
    group_by(rsid) %>% slice_head(n=1) %>% ungroup()
  #There were 9 cases with duplicates- here just picking one arbitrarily at this point
  still.an.issue <- solved %>% group_by(hg38) %>% summarize("count"=n()) %>% dplyr::filter(count > 1)
  stopifnot(nrow(still.an.issue) == 0)
  #This leaves  around 500 SNPs where the top scorer isn't the closest gene by gencode.
  #In this case, we just want to pick theone that is closest based on our coordinates
  #errmm still some duplicates wth

}
  #This takes a matrix of snps that are duplicated, and just returns 
  #(dplyr::filter(snp_gene_map.singular, hg38 %in% multi.snps$hg38),tss_coords

#' List of SNPs that have inconclusive data from OpenTargets (either an NA assignment or tie in the closest genes with no distance measurements)
#' if its tie, I calcualte the distance measures with my reference
#' if its an NA (most likely a bug from extracting from the parquet files), i just map to to the nearest gene with my file
#' @param multi.snp.dat 
#' @param tss_coords 
#'
#' @return
#' @export
#'
#' @examples
tiebreakers_by_distance <- function(multi.snp.dat,tss_coords)
{
  if(!("chr" %in% colnames(multi.snp.dat)))
  {
    multi.snp.dat %<>% separate(hg38,into = c("chr","pos"),sep=":",remove = FALSE)
  }
  tss.coords <- fread(tss_coords) %>% set_colnames(c("chr","start","end","gene", "score","strand")) %>% 
    mutate("hg38" = paste0(gsub(x=chr,pattern = "chr",replacement = ""),":",end)) %>% dplyr::select(-score)
  ret.mat <- NULL
  for(snp in unique(multi.snp.dat$hg38))
  {
    
    cd <- multi.snp.dat %>% dplyr::filter(hg38 %in% snp)
    stopifnot(length(unique(cd$overall_score))==1)
    chr_ <- paste0("chr", unique(cd$chr))
    pos_ = as.numeric(unique(cd$pos))
    top_hits <-(tss.coords %>% dplyr::filter(chr==chr_) %>% mutate("distance"=abs(end-pos_)) %>% arrange(distance))
    matched.genes <- dplyr::filter(top_hits, gene %in% cd$gene_id)
    if(nrow(matched.genes) == 0)
    {
      if(nrow(cd) == 1 & is.na(cd$gene_id))
      {
        message("No gene is NA, using distance")
        ret.hit <- cd[1,]
        ret.hit$d <- top_hits$distance[1]
        ret.hit$gene_id <- top_hits$gene[1]
      }else{
        message("No gene matches, just going with the top 1....")
        ret.hit <- cd[1,]
      }

    }else
    {
      ret.hit <- cd %>% dplyr::filter(gene_id == matched.genes$gene[1])
      ret.hit$d <- matched.genes$distance[1]
    }
    ret.mat <- rbind(ret.mat, ret.hit)
  }
  ret.mat %>% dplyr::select(-chr,-pos)
}
  
  
mapGeneMinDist <- function(df,tss.coords)
  {
  if(!("snp_chr"%in% colnames(df)))
  {
    df$snp_chr <- paste0("chr", sapply(str_split(df$hg38,pattern=":",n=2), function(x) x[1]))
  }
    same.chr <- tss.coords %>% dplyr::select(-hg38) %>% dplyr::filter(chr %in% unique(df$snp_chr)) %>% dplyr::filter(gene %in% df$gene_id) 
    #Make sure the order is the same.....
    if(nrow(same.chr) ==0) #none of the genes offered match... could happen
    {
      message("No gene matches, just going with the top 1....")
      return(df[1,])
    }
    same.chr <- same.chr %>% group_by(gene) %>% mutate("max_end"=max(end)) %>% dplyr::select(chr,gene,max_end) %>% distinct() %>%
      rename(gene="gene_id")
    final.goal <- (left_join(same.chr,df,by="gene_id") %>% mutate("dist"=abs(max_end-snp_end)) %>% arrange(dist))
    ret <- df[1,]; ret$gene_id=final.goal$gene_id[1]
    ret
  }

 

#write_tabular_reports(snp_scores,snp_ids, prioritized_genes, joined.dat,joined.singular, bg_genes, odir)
write_tabular_reports <- function(snp_scores,snp_ids, prioritized_genes, prioritized_snps, joined.dat,joined.singular, bg_genes, odir )
{
  #write out things
  for(i in 1:ncol(snp_scores))
  {
    factor <- gsub(colnames(snp_scores)[i],pattern = "U", replacement = "F")
    test_set <- data.frame("genes"= prioritized_genes[[i]])
    
    combined <- data.frame("V4" = snp_ids, "U"=snp_scores[,i]) %>% left_join(., joined.dat, by="V4") %>%
      arrange(-abs(U)) %>% mutate("rank"=row_number(),"hg37"=gsub(paste0(V1,":",V3),pattern = "chr", replacement = "")) %>%
      dplyr::select(V4,hg38,hg37,gene_id,U,rank) %>% set_colnames(c("RSID","hg38", "hg37","gene", "U", "rank")) %>% rowwise() %>%
      mutate("top_gene"=ifelse(gene %in% prioritized_genes[[i]], TRUE,FALSE)) %>% ungroup() %>%
      mutate("top_snp" =ifelse(RSID %in% prioritized_snps[[i]], TRUE,FALSE))
    
    
    combined.singular <- data.frame("rsid" = snp_ids, "U"=snp_scores[,i]) %>% left_join(., joined.singular, by="rsid") %>%
      arrange(-abs(U)) %>% mutate("rank"=row_number()) %>%
      dplyr::select(rsid,hg38,gene_id,U,rank) %>% set_colnames(c("RSID","hg38", "gene", "U", "rank")) %>% rowwise() %>%
      mutate("top_gene"=ifelse(gene %in% prioritized_genes[[i]], TRUE,FALSE)) %>% ungroup() %>%
      mutate("top_snp" =ifelse(RSID %in% prioritized_snps[[i]], TRUE,FALSE))
    
    bg_set <- data.frame("genes"= bg_genes[[i]])
    count.by.gene.in.set <- data.frame(table(prioritized_genes[[i]])) %>% set_colnames(c("Gene", "Count_in_set")) %>%
      arrange(-Count_in_set)
    #Visualize this
    hist <- ggplot(count.by.gene.in.set, aes(x=Count_in_set)) + geom_histogram(bins=30) + xlab("# of times a gene appears\nin set") + theme_classic(13)
    bar <- ggplot((count.by.gene.in.set %>% dplyr::arrange(-Count_in_set))[1:50,], aes(x=Count_in_set,y=reorder(Gene, Count_in_set))) +
      geom_bar(stat="identity") + xlab("Gene indcidence in set") + theme_classic(13) + ylab("Gene name")
    cowplot::plot_grid(plotlist = list(hist,bar), nrow = 1)
    
    ggsave(paste0(odir, factor,"_gene_freq.png"),height = 10,width=7)
    
    # Write output files, one for each factor...
    write.csv(count.by.gene.in.set, paste0(odir, factor,"_genes_in_factor_count.csv"), row.names = FALSE, quote = FALSE)
    write.table(test_set, paste0(odir, factor,"_genes_in_factor.txt"), row.names = FALSE, quote = FALSE,col.names=FALSE)
    write.table(bg_set, paste0(odir,factor,"_background_genes.txt"), row.names = FALSE, quote = FALSE,col.names=FALSE)
    write.csv(combined, paste0(odir,factor,"_snp_gene_rank.csv"), row.names = FALSE, quote = FALSE)
    write.csv(combined.singular, paste0(odir,factor,"_snp_gene_rank_singular.csv"), row.names = FALSE, quote = FALSE)
  }
}

#This is for snps which wec ouldn't map to hg38
recoverUnmappedSNPs <- function(missing.snps)
{
  #Recover those SNPs with no mapping
  hg37.genes <- fread("/data/abattle4/aomdahl1/reference_data/gencode/gencode.v19.genes.TSS.bed") %>%
    set_colnames(c("chr","start","end","gene","score","strand"))
  hg37.snps <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz") %>% dplyr::filter(rsid %in% unlist(missing.snps))
  #TODO- store the RSID as well as hg38
  added.snps <- NULL
  for(i in 1:nrow(hg37.snps))
  {
    curr.chr <- paste0("chr",hg37.snps[i,]$chrom)
    top.gene <- (dplyr::filter(hg37.genes, chr ==curr.chr ) %>% mutate("dist"=abs(end -hg37.snps[i,]$pos )) %>% arrange(dist))[1,]
    added.snps <- rbind(added.snps, 
                        data.frame("hg38"=NA, "rsid"=hg37.snps[i,]$rsid, "ref"=hg37.snps[i,]$ref, 
                                   "alt"=hg37.snps[i,]$alt, "overall_score"=NA, "gene_id"=top.gene$gene, "distance"=top.gene$dist))
  }
  still.missing.snps <- which(!(missing.snps %in% added.snps$rsid))
  stopifnot(length(still.missing.snps) == 0) #Every SNP got mapped to a gene!
  added.snps
}




#opt$snp_gene_map,opt$snp_id_map, opt$snp_scores, opt$method, opt$output
#snp_gene_map_file <- opt$snp_gene_map; snp_id_map <- opt$snp_id_map; snp_scores_file <- opt$snp_scores; method <- opt$method; odir <- opt$output
main <- function(snp_gene_map_file,snp_id_map, snp_scores_file, method, odir) {
  #FIRST STEPS= MAP SNPS TO GENES
  # Read input files
  message("Loading in reference data (includes openTargets distances)")
  message(snp_gene_map_file)
  snp_gene_map <- data.frame(loadSNPtoGeneMap(snp_gene_map_file))
  
  #Upgraded version of joined singular.. a little more straightforward
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






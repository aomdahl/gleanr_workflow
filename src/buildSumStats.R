###
#Ashton Omdahl
#Script to take matrix of learned loadings and convert these into LDSC summary statistics file format.
#Input required is the (projected) loadings, the number of SNPs in the samples, and the list of SNPs you are using. (recommend hapmap 3)
###

pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr, optparse, magrittr)
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/formatting_functs.R")
option_list <- list(
make_option("--projected_loadings", type = 'character', help = "projected_loadings file. First column is ecpected to be snp ids,  in form chr:pos or RSID"),
make_option("--factors", type = 'character', help = "factor matrix; needed for good estimate of sample sizes."),
make_option("--output", type = "character", help = "Path to output location."),
make_option("--hapmap_list", type = "character", help = "Path SNP list to use, default given.", 
            default = "/data/abattle4/aomdahl1/reference_data/hapmap_chr_ids.txt"),  
make_option("--samp_counts", type = "character", help = "Number of samples across the studies- give as an average.", default = "avg"),
make_option("--samp_file", type = "character", help = "Option to read in a file that contains the sample counts"),
#make_option("--no_file_headers", type = "logical", help = "Should we read in a file header", default = FALSE, action="store_true"),
make_option("--normal_transform", type = "character", help = "Put these on the scale of regular z-scores, so mean is 0, variance is 1", action = "store_true", default = FALSE)
)
#args <- parse_args(OptionParser(option_list=option_list),args=t)
args <- parse_args(OptionParser(option_list=option_list))

#DEBUG:
t <- c("--projected_loadings=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/latent.loadings.txt",
       "--factors=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/latent.factors.txt",
       "--output=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/loading_ss_files/",
       "--samp_file=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.n.tsv")
if (FALSE)
{
  args <- list()
  args$projected_loadings <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results//ukbb_GWAS_h2-0.1_rg-0.9/GLEANER_WL_0.7/no_lambda_correction-K2/projectedMETA_hapmap3_loadings.txt"
  args$samp_file <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results//ukbb_GWAS_h2-0.1_rg-0.9/full_hapmap3_snps.N.tsv"
  args$hapmap_list <- "/data/abattle4/aomdahl1/reference_data/hapmap_chr_ids.txt"
  #Ah, we need the full hapmap list because it has the ref/alt alleles. This is what we are
  args$hapmap_list <- "/data/abattle4/aomdahl1/reference_data/hapmap_chr_ids.txt"
  args$samp_counts <- -1
  args$factors <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results//ukbb_GWAS_h2-0.1_rg-0.9/GLEANER_WL_0.7/no_lambda_correction-K2/latent.factors.SUB.txt"
  args$output <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results//ukbb_GWAS_h2-0.1_rg-0.9/GLEANER_WL_0.7/loading_ss_files_META/"
  factor_names <- c("test1", "test2")
}


#1) get in the projected loadings table
projected_loadings <- fread(args$projected_loadings,header = TRUE)
if(names(projected_loadings)[1] == "chr" && names(projected_loadings)[2] == "pos") #we use the chr/pos version
{
  projected_loadings[,1] <- paste0(projected_loadings$chr, ":", projected_loadings$pos)
  projected_loadings <- projected_loadings %>% select(-pos) %>% rename("id" = chr)
  #map_names <- c("id", "RSID")
  mapping <- fread(args$hapmap_list)  %>% mutate(id = paste0(`CHROM`, ":", POS))
  dat <- left_join(projected_loadings, mapping, by = "id") %>% arrange(id) %>% 
    select(-id) %>% rename("id" = ID) %>% arrange(id)
} else if(grepl(pattern = ":",x = projected_loadings[1,1])) { #version based on SNP ID...
  colnames(projected_loadings)[1] = "id"
  mapping <- fread(args$hapmap_list)  %>% mutate(id = paste0(`CHROM`, ":", POS))
  dat <- left_join(projected_loadings, mapping, by = "id") %>% arrange(id) %>%
    select(-id) %>% rename("id" = ID) %>% arrange(id)
}else { #The rsid version..
  update_names <- names(projected_loadings)
  update_names[1] <- "id"
  names(projected_loadings) <- update_names
  #map_names <- c("chr_pos", "id")
  mapping <- fread(args$hapmap_list)  %>% rename("id" = ID)
  dat <- left_join(projected_loadings, mapping, by = "id") %>% arrange(id)
}
print("projected_loadings inputted and processed")
#factor_names <- scan(args$factor_names, what = "character")
#c1: 1:pos
#c2:factor....

#CHROM  POS ID  REF ALT

#names(mapping) <- map_names


#2) get the snp order

print("Ids lined up!")

#2) get the counts the right way
if(args$samp_counts == "avg" | args$samp_counts == "weighted")
{
  #a few options for doing this- take the average across all studies, the sum, etc.
  #For now, I'm not sure what to do. but the weighting doesn't seem like its the right way.
  n <- fread(args$samp_file) %>% drop_na() 
  nn <- standardizeColNames(colnames(n)); nn[1] <- "ids";
  n <- n %>% set_colnames(nn) %>% filter(ids %in% dat$id) %>% arrange(ids)
  if(args$samp_counts == "avg")
  {
	  message("Sample size by mean")
	  var_counts <- rowMeans(n[,-1])
  }
  if(args$samp_counts == "weighted")
  {
	  #12/02- better heuristic for weighting...
	  #Need to make sure the names align. That's the point...
	  factors <- fread(args$factors) %>% mutate("rownames" = factor(rownames, level=nn)) %>% arrange(rownames)
	  stopifnot(all(nn[-1] == factors$rownames))
	  #order rows 
	  weights <- apply(factors[,-1],2, function(x) x^2/sum(x^2))
	  weighted.n <- as.matrix(n[,-1]) %*% as.matrix(weights)
	  #var_counts <- rowMeans(n[,-1])
	  if(nrow(weighted.n) != nrow(dat))
	  {
	    message("More SNPs to evaluate than we've sampled SNPs for.")
	    message("Current standard- reduce L_proj to the smaller size.")
	    ndat <- filter(dat,id %in% n$ids)
	    if(nrow(ndat) > 900000)
	    {
	      dat <- ndat
	    }
  	  }
  }
  }  else {
  var_counts = rep(as.numeric(args$samp_counts), nrow(dat))
}
#3) build a sumstats file for each factor ''
offset <- 1

if( !grepl(":", dat[1,1]) & !grepl("rs", dat[1,1]))
{
  message("DATA NOT CORRRECTLY FORMATTED, ENDING NOW.")
  quit()
}

if( !grepl(":", projected_loadings[1,1]) & !grepl("rs", projected_loadings[1,1]))
{
  message("DATA NOT CORRRECTLY FORMATTED, ENDING NOW.")
  quit()
}


for(i in 2:ncol(projected_loadings)) #from 2 since first column is IDs...
{
  #SNP   A1  A2  N   Z #a1 is effect allel
  if(!(names(projected_loadings)[i] %in% c("id", "CHROM", "POS", "REF", "ALT"))) #make sure its the factr columns
  {
    
    #this is a bit.... concerning. Why the error????
    zs <- as.numeric(unlist(dat[,..i]))
  print(head(zs))
    #Briefly report on the distribution of the data.x
    rand.samp <- sample(zs, size = 10000)
    #message("Data has mean ",   round(mean(rand.samp, na.rm = TRUE), digits = 5), " and variance ", round(var(rand.samp, na.rm = TRUE), digits = 5))
    #dev = ks.test(rand.samp, y="pnorm")$p.value
    #message("By the KS test of deviation from standard normal, pvalue is: ", dev)
    if(args$normal_transform)
    {
      zs <- scale(zs) #This is not an inverse rank normal transform. That may be better.
    }
    if(args$samp_counts == "weighted")
    {
	    var_counts = weighted.n[,(i-offset)]
    }
    out <- data.frame("SNP" = dat$id, "A1" = dat$ALT, "A2" = dat$REF, "N" = var_counts, "Z" = zs) #TODO: check the ref alt stuff make sure on point.
    names(out) <- c("SNP", "A1","A2", "N","Z")
    #write_tsv(out, paste0(args$output,"/F",i-offset, ".sumstats.gz"))
    fwrite(x = out, file= paste0(args$output,"/F",i-offset, ".sumstats.gz"), compress = "gzip", sep = "\t")
    message(paste0("written out ", file.path(args$output), "/F", i-offset, ".sumstats.gz"))
  } else {offset <- offset + 1}

}

#4) run the ldsc for each one!

###
#Ashton Omdahl
#Script to take matrix of learned loadings and convert these into LDSC summary statistics file format.
#Input required is the (projected) loadings, the number of SNPs in the samples, and the list of SNPs you are using. (recommend hapmap 3)
###

matchOrder <- function(reference, alternate)
{
  #clean up stupid artefacts:
  reference <- sapply(reference, function(a) gsub(pattern="^X_",replacement = "_", x=a))
  quick.check <- sapply(1:length(reference), function(i) grepl(pattern =reference[i], x=alternate[i]))
  if(all(quick.check))
  {
    return(1:length(alternate))
  }
  for(r in reference)
  {
    ret <-c()
    lookups <- which(sapply(alternate, function(x) grepl(pattern = r,x = x)))
    if(length(lookups) > 1 | length(lookups) < 1 )
    {
      message("Multiple matches or too few matches- evaluate")
      message("Adding a - before and after to see if that works.")
      lookups <- which(sapply(alternate, function(x) grepl(pattern = paste0("-",r, "-"),x = x)))
      
      if(length(lookups) > 1 | length(lookups) < 1 )
      {
        message("Still no success")
        print(r)
        print(lookups)
        break()
      }
    }else
    {
      ret <- c(ret,lookups)
    }

  }
  message("completed with", r)
  ret
}

pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr, optparse, magrittr)
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/formatting_functs.R")
option_list <- list(
make_option("--loadings", type = 'character', help = "loadings file. First column is ecpected to be snp ids,  in form chr:pos or RSID"),
make_option("--factors", type = 'character', help = "factor matrix; needed for good estimate of sample sizes."),
make_option("--output", type = "character", help = "Path to output location."),
make_option("--samp_counts", type = "character", help = "Number of samples across the studies- give as an average.", default = "weighted"),
make_option("--samp_se", type = "character", help = "SE OF SNPS ON INPUT."),
make_option("--samp_file", type = "character", help = "Option to read in a file that contains the sample counts"),
make_option("--no_zscaling", type = "logical",default=FALSE, action="store_true", help = "Specify this if you don't want to scale by SEs (for SVD or flash)"),
make_option("--snp_reflist", type = "character", help = "Path SNP list to use, default given.", 
            default = "/data/abattle4/aomdahl1/reference_data/hapmap_chr_ids.txt")
#make_option("--no_file_headers", type = "logical", help = "Should we read in a file header", default = FALSE, action="store_true"),
)
t <- c("--loadings=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/latent.loadings.txt",
       "--factors=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/latent.factors.txt",
       "--output=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/loading_ss_files/",
       "--samp_file=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.n.tsv",
       "--samp_se=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.1.se.tsv",
       "--snp_reflist=/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz")
#args <- parse_args(OptionParser(option_list=option_list),args=t)
args <- parse_args(OptionParser(option_list=option_list))

#DEBUG:

#1) get in the projected loadings table
loadings <- fread(args$loadings,header = TRUE)
update_names <- names(loadings)
update_names[1] <- "id"
names(loadings) <- update_names
loadings <- loadings %>% arrange(id)
print("loadings inputted and processed")

#2) get the snp order
print("Ids lined up!")

#2) get the counts the right way
  #a few options for doing this- take the average across all studies, the sum, etc.
  #For now, I'm not sure what to do. but the weighting doesn't seem like its the right way.
  n <- fread(args$samp_file) %>% filter(ids %in% loadings$id)  
  #Set NAs to 0, that means that study isn't contributing at that SNP most likely...
  nn <- standardizeColNames(colnames(n)); nn[1] <- "ids";
  n <- n %>% set_colnames(nn) %>% filter(ids %in% loadings$id) %>% arrange(ids)
  
  ###Now get the se
  se <- fread(args$samp_se) %>% filter(ids %in% loadings$id)  
  sen <- standardizeColNames(colnames(se)); sen[1] <- "ids";
  se <-se %>% set_colnames(sen) %>% filter(ids %in% loadings$id) %>% arrange(ids)
  
  stopifnot(all(se$ids == n$ids))
  if(args$samp_counts == "avg")
  {
	  message("Sample size by mean")
	  var_counts <- rowMeans(n[,-1])
  }
  if(args$samp_counts == "weighted")
  {
	  #12/02- better heuristic for weighting...
	  #Need to make sure the names align. That's the point...
	  factors <- fread(args$factors) 
	  
	  #Sample sizes
	  new.order <- matchOrder(factors$Study, nn[-1]) + 1
	  n <- n[,..new.order] %>%  apply(., 2, as.numeric)
	  n[is.na(n)] <- 0 #set those eto 0
	  
	  #SE
	  new.order <- matchOrder(factors$Study, sen[-1]) + 1
	  se <- se[,..new.order] %>%  apply(., 2, as.numeric)
	  w <- 1/se
	  w[is.na(w)] <- 0 #set those eto 0
	 #%>% mutate("rownames" = factor(rownames, level=nn)) %>% arrange(rownames)
	  #order rows 
	  weights <- apply(factors[,-1],2, function(x) x^2/sum(x^2))
	  weighted.n <- as.matrix(n) %*% as.matrix(weights)
	  weighted.se <- as.matrix(w) %*% as.matrix(weights)
	  #var_counts <- rowMeans(n[,-1])
	  if(nrow(weighted.n) != nrow(loadings))
	  {
	    message("More SNPs to evaluate than we've sampled SNPs for.")
	    message("Current standard- reduce L_proj to the smaller size.")
	    ndat <- filter(loadings,id %in% n$ids)
	    if(nrow(ndat) > 900000)
	    {
	      loadings <- ndat
	    }
  	 }
  } 

#3) build a sumstats file for each factor ''

if( !grepl(":", loadings[1,1]) & !grepl("rs", loadings[1,1]))
{
  message("DATA NOT CORRRECTLY FORMATTED, ENDING NOW.")
  quit()
}

dat <- as.matrix(loadings[,-1])

#Process the SNP ref list:
snp.dat <- fread(args$snp_reflist)
colnames(snp.dat) <- toupper(colnames(snp.dat))
if("RSID" %in% colnames(snp.dat))
{
  snp.dat <- snp.dat %>% rename("id"=RSID)
}
snp.dat <- filter(snp.dat, id %in% loadings$id) %>% arrange(id)
stopifnot(all(snp.dat$id == loadings$id))
message("Weighting by the inverted SE")
for(i in 1:ncol(dat)) #from 2 since first column is IDs...
{
  #SNP   A1  A2  N   Z #a1 is effect allele....
    zse <- dat[,i] * weighted.se[,i]
    if(args$no_zscaling)
    {
      zse <- dat[,i]
    }
    zs.n <- dat[,i] * weighted.n[,i]
    #Briefly report on the distribution of the data.x
    #rand.samp <- sample(zse, size = 10000)
    #message("Data has mean ",   round(mean(rand.samp, na.rm = TRUE), digits = 5), " and variance ", round(var(rand.samp, na.rm = TRUE), digits = 5))


    if(args$samp_counts == "weighted")
    {
	    var_counts = weighted.n[,i]
    }
    message("We assume the ALT is the effect allele, ensure this is true for your reference and dataset")
    out <- data.frame("SNP" = snp.dat$id, "A1" = snp.dat$ALT, "A2" = snp.dat$REF, "N" = var_counts, "Z" = zse) #TODO: check the ref alt stuff make sure on point.
    names(out) <- c("SNP", "A1","A2", "N","Z")
    fwrite(x = out, file= paste0(args$output,"/F",i, ".sumstats.gz"), compress = "gzip", sep = "\t")
    
    #alternative N weighted version
    #Don't make this anymore
    #out <- data.frame("SNP" = snp.dat$id, "A1" = snp.dat$ALT, "A2" = snp.dat$REF, "N" = var_counts, "Z" = zs.n) #TODO: check the ref alt stuff make sure on point.
    #names(out) <- c("SNP", "A1","A2", "N","Z")
    #fwrite(x = out, file= paste0(args$output,"/F",i, "N_weighted.sumstats.gz"), compress = "gzip", sep = "\t")
    
    
    message(paste0("written out ", file.path(args$output), "/F", i, ".sumstats.gz"))

  } 




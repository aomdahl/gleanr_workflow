args = commandArgs(TRUE)
#commands are 
#1) seed
#2) # simulations to run/number of iterations to run
#3) the  directory name you want to write in (a folder that should exist.)
#4) The diriectory path to write in (sans #5 above)
#5) Heritability of trait 1
#6) Heritability of trait 2
#7) Noise on the 3 phenotyppes (default is 1,1,1)
#8) Overlap list (default is c(0,1000,5000,10000,15000,20000,25000,30000))
#9) The overall sample size
#Rscript simulation_scripts/gwas_centric_vary_overlap_degree_template.R 5 10 large_studies_with_covar  /scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/ 0.12 0.05
#args <- c("5", "10", "large_studies_with_covar",
#          "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/",
#          "0.12", "0.05")
#args <- c("5", "10", "large_studies_with_covar", "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/", "0", "0", "1")
#debug
#args <- c("23", "10", "$ODIR",  "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/null_matrix_spearman/", "0", "0", "1", "0,10000,15000,20000,23000,26000,30000", "30000")

set.seed(args[1])
seedi = args[1]
herit_a=as.numeric(args[5])
message("input arguments are...")
print(args)
if(herit_a == -1)
{
  message("going to make all the traits null")
  herit_b=0
  herit_a=0
}else
{
  herit_b=as.numeric(args[6])
}

message("Seed set to ", args[1])
message("this version of the script simply saves out all the GWAs summary stats, doesn't actually do any analysis.")
message("herit a: ", herit_a)
message("herit b: ", herit_b)

message("Number of simulations to run:", args[2])
n.iter <- as.numeric(args[2])
dir_tag <- args[3]
odir = args[4]
#processing the phenotype thing
esd <- as.numeric(stringr::str_split(args[7],pattern = ",")[[1]])
message("SD of noise to add to sims: ", esd[1], ", ", esd[2], ", ", esd[3] )

#Processing overlap list
overlap.list <- c(0,1000,5000,10000,15000,20000,25000,30000)
if(length(args) >= 8)
{
  message("Processing overlap list")
  overlap.list <- as.numeric(stringr::str_split(args[8],pattern = ",")[[1]])
  message("Custom overlap list: " ,paste(overlap.list, collapse = ","))
}
#odir="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/gwas_sim_style/non_null_matrix/"
pacman::p_load(stats, data.table, magrittr, dplyr, tidyr, ggplot2, stringr, fossil)
source("src/helper_functions_plieo.R")
#source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source('../../../../src/plot_functions.R')
#load("/scratch16/abattle4/ashton/snp_networks/scratch/gwasMF_pleiotropy/.RData")

#Generate genotypes for 3 populations...
## These are based on MAFs from the thousand genomes, subsetted to a pruned set used in real analysis
message("Loading MAF data")
thou.g.mafs <- fread("/scratch16/abattle4/ashton/prs_dev/1000genomes_refLD/plink.frq")
#thou.g.pruned <- fread("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/500kb.0.04r2.prune.in", header = F)
#thou.g.pruned <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark/ukbb_benchmark.250kb.0.2r2.prune.in", header = F)
thou.g.pruned <- fread("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2/p1e-5/ukbb_benchmark_2.250kb.0.2r2.prune.in", header = F)
thou.g.mafs.sub <- thou.g.mafs %>% filter(SNP %in% thou.g.pruned$V1)
thou.g.mafs.use <- thou.g.mafs.sub
rm(thou.g.mafs)

#First- create 3 LARGE cohorts of 100,000 individuals each.
if(length(args) >=9)
{
	n.fixed = as.numeric(args[9])
	N_tot = n.fixed * 3 + 1000
}else
{
N_tot <- 100000
n.fixed = 30000
}
message("Populations will have ", N_tot, " for ", n.fixed, " samples per study")

#1
genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
pop1.geno <- t(do.call("cbind", genotypes.all))
#2
genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
pop2.geno <- t(do.call("cbind", genotypes.all))

#3
genotypes.all <- lapply(1:N_tot, function(x) genGenotypes(nrow(thou.g.mafs.use), thou.g.mafs.use$MAF))
pop3.geno <- t(do.call("cbind", genotypes.all))
rm(genotypes.all)


#Simulation settings
rep.dat <- list()
message("Starting sim iters...")
n.snps <- 1000
snp.weights <- list()
snp.weights[[1]] <- c(rnorm(200,sd=sqrt(herit_a/200)), rep(0,800))
snp.weights[[2]] <- c(rep(0,200), rnorm(300,sd=sqrt(herit_b/300)),rep(0,500))
snp.weights[[3]] <- c(rep(0,10), rep(0,900),rep(0.07,90))
stopifnot(length(snp.weights[[3]]) == n.snps)
stopifnot(length(snp.weights[[2]]) == n.snps)
stopifnot(length(snp.weights[[1]]) == n.snps)
n.groups=3
n.phenos = 3
var.type="1+2=3"
if(n.phenos == 5)
{
  var.type="1+2=3,4~5"
}
if(length(esd) == 1)
{
  message("going for BMi-like, null")
  var.type="bmi_like"
}
message("Starting the loop")
#Get the phenotypes
for(rep in 1:n.iter)
{
    start = Sys.time() 
    message("on iteration ", rep, " of ", n.iter)
    #Simulate new phenotypes
  
    pop1.null <- list()#genotypes.tab
 if(herit_a==0 & herit_b==0)
 	{   
		message("In the null builder")
   #Correction- esd was going through wrong....
    y1 <- phenotypeBuilder(NULL, pop1.geno[,1:n.snps], esd, type = var.type)
    y2 <- phenotypeBuilder(NULL, pop2.geno[,1:n.snps], esd, type = var.type)
    y3 <- phenotypeBuilder(NULL, pop3.geno[,1:n.snps], esd, type = var.type)
	}else
	{
    y1 <- phenotypeBuilder(snp.weights, pop1.geno[,1:n.snps], esd, type = var.type)
    y2 <- phenotypeBuilder(snp.weights, pop2.geno[,1:n.snps], esd, type = var.type)
    y3 <- phenotypeBuilder(snp.weights, pop3.geno[,1:n.snps], esd, type = var.type)	

	}
    pop2.null <- list() #genotypes.tab.2
    pop3.null <- list()
  
    first <- list(); second <- list(); third <- list();
    overlap.at.level <- list()
    full_ret = TRUE
    #NOTE: index is starting at 1, not corresponding to overlap list value. This changes elsewhere, so beware.
    for(i in 1:length(overlap.list))
    {
      print(i)
      l <- overlap.list[[i]]
      #calcCovar <- function(n.snps, n.overlap, n.groups,n.pheno=3, y.list, N.tot = 10000, fixed.count = 5000)
      overlap.at.level[[as.character(l)]] <- calcCovar(n.snps, l, n.groups, n.pheno = n.phenos,y.list = list(y1,y2,y3), N.tot=n.fixed)

      first[[as.character(l)]] <- relatedOverlapBETTER(n.snps, l, 3, y1, N.tot = N_tot,
                                             genotypes=pop1.geno[,1:n.snps], fixed.count = n.fixed,ret.all.stats=full_ret)

     second[[as.character(l)]] <- relatedOverlapBETTER(n.snps, l, 3, y2, N.tot = N_tot,
                                             genotypes=pop2.geno[,1:n.snps], fixed.count = n.fixed,ret.all.stats=full_ret)

      third[[as.character(l)]] <- relatedOverlapBETTER(n.snps, l, 3, y3, N.tot = N_tot,
                                             genotypes=pop3.geno[,1:n.snps], fixed.count = n.fixed,ret.all.stats=full_ret)


		
     
    }
    
    #Save the data with different combinations of overlap for each trait
    all.options <- expand.grid(overlap.list, overlap.list, overlap.list)
    #simplify this- we just need one level of overlap per setting, 
    #i.e. we don't need c(10,0,0), c(0,10,0), and (0,0,10), these should be the same
    all.ids <- sapply(1:nrow(all.options), function(i)  paste(sort(unlist(all.options[i,])), collapse="-"))
    sub.options <- all.options[!duplicated(all.ids),]
  
      for(opt_i in 1:nrow(sub.options))
      {
        opt=sub.options[opt_i,]
        message("Current overlap setting is: ", opt)
      #all the same
      f <- as.character(opt[1]);s <- as.character(opt[2]);t <- as.character(opt[3])
      #Get the cohort overlap things....
      f_o <- overlap.at.level[[f]]$sample.cor[[1]]; s_o <- overlap.at.level[[s]]$sample.cor[[2]]; t_o <- overlap.at.level[[t]]$sample.cor[[3]];
      #get the SE
      f_se <- overlap.at.level[[f]]$sample.se[[1]]; s_se <- overlap.at.level[[s]]$sample.se[[2]]; t_se <- overlap.at.level[[t]]$sample.se[[3]];
      #each object contains 2 lists, each with 3 matrices
      covar.overall <- as.matrix(Matrix::bdiag(list(f_o, s_o, t_o)))
      se.overall <- as.matrix(Matrix::bdiag(list(f_se, s_se, t_se)))
      #Write each out
      write.table(data.frame(round(cbind(first[[f]]$se, second[[s]]$se, third[[t]]$se),digits = 6)) %>% set_colnames(c("A1", "B1", "C1","A2", "B2", "C2","A3", "B3", "C3")),  
                  file = paste0(odir, dir_tag, "/seed", seedi,".overlap-",f,"-",s,"-",t, ".replicate",rep,".SE.csv"),
                  quote = FALSE,sep = ",", row.names = FALSE )
      write.table(data.frame(round(cbind(first[[f]]$beta, second[[s]]$beta, third[[t]]$beta),digits = 6)) %>% set_colnames(c("A1", "B1", "C1","A2", "B2", "C2","A3", "B3", "C3")),  
                  file = paste0(odir, dir_tag, "/seed", seedi,".overlap-",f,"-",s,"-",t, ".replicate",rep,".BETA.csv"),
                  quote = FALSE,sep = ",", row.names = FALSE )
      #Covar details
      write.table(data.frame(covar.overall) %>% set_colnames(c("A1", "B1", "C1","A2", "B2", "C2","A3", "B3", "C3")),  
                  file = paste0(odir, dir_tag, "/seed", seedi,".overlap-",f,"-",s,"-",t, ".replicate",rep,".COV.csv"),
                  quote = FALSE,sep = ",", row.names = FALSE )
      write.table(data.frame(se.overall) %>% set_colnames(c("A1", "B1", "C1","A2", "B2", "C2","A3", "B3", "C3")),  
                  file = paste0(odir, dir_tag, "/seed", seedi,".overlap-",f,"-",s,"-",t, ".replicate",rep,".COV_SE.csv"),
                  quote = FALSE,sep = ",", row.names = FALSE )
      
      }

}


source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/selectivePressureEst.R")
pacman::p_load(magrittr, dplyr, ggplot2, data.table,optparse)

#' Calculate the selection scores
#'
#' This function calculates the selection scores based on provided U matrices and processed data.
#'
#' @param u.mat The U matrix.
#' @param u.processed U data, including MAF and other relevant scores. Lazy way of doing this
#' @param permute Number of permutations to run, default is 0
#' @return A dataframe with calculated selection scores.
#' @export
getSelectionScores <- function(u.mat,u.processed, permute = 0)
{
  u.selection.full.permutes <- list()
  u.pvals.s_sigma <-  matrix(NA, ncol = 2, nrow = ncol(u.mat))
if(permute > 0)  u.pvals.s_sigma <-  matrix(NA, ncol =4, nrow = ncol(u.mat))
  
  for(u in 1:ncol(u.mat))
  {
    u_i = u+1
    #Select out the non-zero effects
    u_effects <- u.processed[,u_i][u.processed[,u_i] != 0]
    #And the corresponding minor allele frequencies
    mafs <-  u.processed$maf[u.processed[,u_i] != 0]
    u.optim <- estimateS_SigmaG(u_effects, mafs)
    if(permute > 0)
    {
      u.selection.full.permutes[[u]] <-permuteNullDistS(2,10, u_effects, mafs) 
      u.pvals.s_sigma[u,] <- c(u.optim$par[1], u.optim$par[2], sum(u.selection.full.permutes[[u]][,1] <= u.optim$par[1]),
                               sum(u.selection.full.permutes[[u]][,2] >= u.optim$par[2]))
    }else
    {
      u.pvals.s_sigma[u,] <- c(u.optim$par[1], u.optim$par[2])
    }

  }
  out.col.names <- c("S_hat", "sigma_g_hat")
  if(permute > 0) out.col.names <- c("S_hat", "sigma_g_hat", "p_s", "p_sigma_g")
  df.s_estimates <- data.frame(u.pvals.s_sigma) %>% set_colnames(out.col.names) %>% 
    mutate("Factor" = paste0("U", 1:ncol(u.mat))) %>% print()
  
  #Also do factor sparsity,
  df.s_estimates$sparsity <- apply(u.mat, 2, function(x) sum(x==0)/length(x))
  return(df.s_estimates)
}

#' Load Pan-UKBB MAF reference
#'
#' This function loads the Pan-UKBB MAF reference based on provided SNP IDs and file path.
#'
#' @param snp.ids SNP IDs to subset to
#' @param maf.path Path to MAF reference file.
#' @return A dataframe containing filtered Pan-UKBB MAF reference.
#' @export
loadMAFRef <- function(snp.ids, maf_id="af_EUR", maf.path="/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz")
{
  full.snps <- fread(maf.path)
  full.snps <- filter(full.snps,rsid %in% snp.ids) #There are 2 redundant SNPs it looks like.
  i=which(colnames(full.snps) ==maf_id)
  print(i)
  print(full.snps[,..i])
  full.snps$maf <- sapply(unlist(full.snps[,..i]), function(x) ifelse(x > 0.5, 1-x, x)) + 1e-16
  #adding a constant so none are zero, which kills the function
  full.snps
}

#' Join U Data with MAF, aF and rsid
#'
#' This function joins U data with SNP reference data.
#'
#' @param U The U matrix.
#' @param snp.ids SNP IDs.
#' @param snp.ref MAF reference data
#' @return Joined dataframe with additional columns
#' @export
joinedUDat <- function(U, snp.ids, snp.ref)
{
  df.u <- data.frame("snp" = snp.ids, U) %>%  set_colnames(c('rsid', paste0("U",1:ncol(U))))
  df.u <- left_join(df.u, snp.ref, by = "rsid") %>% filter(!is.na(maf))
  
  #There are a no NAs, good job
  if(sum(is.na(df.u$af_EUR)) > 0 |   sum(is.na(df.u$maf)) > 0)
  {
    message( sum(is.na(df.u$maf)), "  of the SNPs don't have mafs, be warned.")
  }
  df.u

}

joinedVDat <- function(vmat, trait.names)
{
  data.frame("trait" =trait.names, vmat) %>% set_colnames(c('trait', paste0("V",1:ncol(vmat))))
}

visualizeSelectionScores <- function(df.s_estimates)
{
 
  f_s <- ggplot(df.s_estimates, aes(x = reorder(Factor,S_hat), y= S_hat, fill = sparsity)) + geom_bar(stat="identity") + 
    coord_flip() + xlab("Factor") + ylab(bquote(hat(S))) + theme_bw() + 
    ggtitle("Selective pressure by factor")
  #S hat wrt sparsity- not really a relationship, with the exception of just one factor.
  s_sparsity <- ggplot(df.s_estimates, aes(x= sparsity,y=S_hat)) + geom_point() + theme_bw() + ggtitle("Relationship between S and factor Sparsity")
  
  
  #S hats wrt sigmas
  s_sigmas <- ggplot(df.s_estimates, aes(x = sigma_g_hat, y= S_hat)) + geom_point() + 
    xlab(bquote(sigma[g])) + ylab(bquote(hat(S))) + theme_bw()  + ggtitle("Relationship between S and sigma_g")
  
  f_sigmas <- ggplot(df.s_estimates, aes(x = reorder(Factor,sigma_g_hat), y= sigma_g_hat)) + geom_bar(stat="identity") + 
    coord_flip() + xlab("Factor") + ylab(bquote(sigma[g])) + theme_bw() + 
    ggtitle(bquote(sigma[g] ~ " by factor"))
  list(f_s, s_sparsity, s_sigmas, f_sigmas)
  
}

plotTraits <- function(factor_list,df.s_estimates,df.v)
{
  plot.list <- list()
  for(q in factor_list)
  {
    q_ <- gsub(x=q, pattern= "U", replacement = "V")
    s <- (df.s_estimates %>% filter(Factor == q))$S_hat
    
    plot.list[[q_]] <- ggplot(df.v %>% filter(abs(!!sym(q_)) > 0), aes(x=reorder(trait, !!sym(q_)), y=!!sym(q_))) + 
      geom_bar(stat="identity") + coord_flip() + theme_bw() + xlab("Traits") + ggtitle(paste0(q, " S=",s))
  }
  plot.list
}

plotMAFvsU <- function(factor_list,df.u)
{
  df.u$maf_term <- df.u$maf*(1-df.u$maf)
  plot.list <- list()
  for(q in factor_list)
  {
    q_i <- which(colnames(df.u) == q)
    plot.list[[q]] <- ggplot(data=df.u, aes(x=maf_term,y=.data[[q]]^2)) + geom_point() + theme_bw() + ylab(bquote(U^2)) + xlab(bquote(f(1-f))) + 
     ggtitle(paste0("Relationship between ", q, " and MAF"))
  }
  plot.list
}

plotAndTestFactorsVsS <- function(V, U, s_dat)
{
  traits.per.factor <- apply(V, 2, function(x) sum(x!=0))
  snps.per.factor <- apply(U, 2, function(x) sum(x!=0))
  message("test of snps per factor vs S")
  #cor.test(s_dat$S_hat, snps.per.factor,method = "spearman")
  print(cor.test(s_dat$S_hat, snps.per.factor,method = "pearson"))
  message("test of traits per factor vs S")
  # cor.test(s_dat$S_hat, traits.per.factor,method = "spearman") 
  print(cor.test(s_dat$S_hat, traits.per.factor,method = "pearson"))
  df.test <- data.frame(df.s_estimates,traits.per.factor,snps.per.factor)
  list(ggplot(df.test, aes(x=traits.per.factor, S_hat)) + geom_point() + theme_bw() + xlab("Traits per factor"),
       ggplot(df.test, aes(x=snps.per.factor, S_hat)) + geom_point() + theme_bw() + xlab("SNPs per factor"))
  
}




saveFigs <- function(ilist, odir, tag)
{
  i=1
  for(image in ilist)
  {
    ggsave(plot = image, filename = paste0(odir, tag, "_", i, ".png"))
    i=i+1
  }
}

# Argument parsing
option_list <- list(
  make_option(c("-f", "--factorization"),
              type = "character", default = NULL,
              help = "Path to an Rdata object of factorization to load in"),
  make_option(c("-a", "--maf_column"),
              type = "character", default = "af_EUR",
              help = "Specify the name of the MAF column to use"),
  make_option(c("-m", "--maf_reference"),
              type = "character", default = "/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz",
              help = "Path to MAF reference file"),
  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = "Output location")
)

opt_parser <- OptionParser(option_list = option_list)
t=c("--factorization=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/_final_dat.RData",
    "--output=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/")
#-a af_EAS -m results//panUKBB_complete_41K/full_maf_dat.txt -o results//panUKBB_complete_41K/selective_pressure_EAS_af/
t=c("--factorization=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K/_final_dat.RData",
    "--maf_column=af_EAS", 
    "--maf_reference=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results//panUKBB_complete_41K/full_maf_dat.txt",
    "--output=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/")
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$factorization) || is.null(opt$output)) {
  stop("Please provide both factorization and output arguments.")
}

# Load necessary data
if (!is.null(opt$factorization)) {
  load(opt$factorization)
}


# Example: 
# results <- your_function(factorization, opt$maf_reference)
message("Loading MAF reference")
full.snps <- loadMAFRef(ret$snp.ids,maf_id=opt$maf_col, 
                        maf.path=opt$maf_reference)
#Plot the MAF
png(paste0(opt$output, "maf_dist_hist.png"))
hist(full.snps$maf, main = "Distribution of MAF in sample", xlab = "MAF", breaks=50,xaxt='n')
axis(side=1, at=seq(from=0, to=0.5, by=0.05))
dev.off()

#Reference data for U
df.u <- joinedUDat(ret$U, ret$snp.ids, full.snps)
df.v <- joinedVDat(ret$V, ret$trait.names)



#Get selection scores:
message("Estimating S scores")
df.s_estimates <- getSelectionScores(ret$U,df.u, permute = 0)
write.table(df.s_estimates, file=paste0(opt$output, "s_scores.tsv"), quote=FALSE, row.names=FALSE)
#Visualize these
all.plots <- visualizeSelectionScores(df.s_estimates)

#Plot the traits in the top 3 factors and bottom 3 factors, as well as relationship between U and MAF
message("Visualizing traits in top and bottom scores")
top.factor.list <- (df.s_estimates %>% arrange(S_hat))[1:3,]$Factor
top.factors <- plotTraits(top.factor.list,df.s_estimates,df.v)
bottom.factor.list <- (df.s_estimates %>% arrange(-S_hat))[1:3,]$Factor
lower.factors <- plotTraits(bottom.factor.list,df.s_estimates,df.v)

#Plot the relationship betwen MAF and Scores for desired traits
top.rel <- plotMAFvsU(top.factor.list,df.u)
lower.rel <- plotMAFvsU(bottom.factor.list,df.u)


#Compare the selective scores to the number of traits loaded on each one

selective.vs.factors <- plotAndTestFactorsVsS(ret$V, ret$U, df.s_estimates)

# Save results
message("Saving out results.")
save(df.s_estimates, all.plots, top.factors, lower.factors, top.rel, lower.rel, file = paste0(opt$output, "/allFigs.RData"))
saveFigs(all.plots, opt$output, "s_scores")
saveFigs(top.factors, opt$output, "topS_traits")
saveFigs(lower.factors, opt$output, "bottomS_traits")
saveFigs(top.rel, opt$output, "topS_maf")
saveFigs(lower.rel, opt$output, "bottomS_maf")
safeFigs(selective.vs.factors, opt$output, "S_vs_factors")
message("Preliminary selective pressure analysis is complete.")

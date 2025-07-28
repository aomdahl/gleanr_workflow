#First the heatmap:
pacman::p_load(data.table, dplyr, magrittr, tidyr)

#Analysis was performed using script `run_scripts/finngen_ukbb_benchmark/v2_expanded_run/covar_influence/finngen_ukbb_benchmark_svd.sh`
#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/R/read_in_tools.R")
#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/R/data_decorrelation.R")
library(gleanr)
#Review this...
alignAndCompare <- function(finngen.svd, ukbb.svd)
{
  v.dat <- prepMatricesForAnalysis(finngen.svd$V,ukbb.svd$V); 
  u.dat <- prepMatricesForAnalysis(finngen.svd$U,ukbb.svd$U)
  stopifnot(!v.dat$swap)
  fg.v <- v.dat$lead; uk.v <- v.dat$second
  fg.u <- u.dat$lead; uk.u <- u.dat$second
  #Get the matchup
  svd.version <- greedyMaxPairedCor(fg.u,fg.v, uk.u,uk.v)
  #align UKBB to match with Finngen
  #uk.u.new = matrixSignsProduct(uk.u[,svd.version$order.pred], svd.version$signs)
  #uk.v.new = matrixSignsProduct(uk.v[,svd.version$order.pred], svd.version$signs)
  #list("adj.v"=cor.test(y=c(uk.v.new), x=c(fg.v)),"adj.u"= cor.test(y=c(uk.u.new), x=c(fg.u)))
  list("adj.v"=svd.version$corr_v,"adj.u"= svd.version$corr_u, x=c(fg.u))
}

getRRMSE <- function(finngen.svd, ukbb.svd)
{
  finngen.x <- as.matrix(finngen.svd$U %*% diag(finngen.svd$d) %*% t(finngen.svd$V))
  ukbb.x <- as.matrix(ukbb.svd$U %*% diag(ukbb.svd$d) %*% t(ukbb.svd$V))
  rrmse(ukbb.x, finngen.x)
}

scaleBoth <- function(x)
{
  x$V <- unitNorms(x$V); x$U <- unitNorms(x$U);
  x
}
######################################################
### Looking directly at correlation:
source("/scratch16/abattle4/ashton/snp_networks/scratch/finngen_v_ukbb/src/matrix_comparison_utils.R")
fg.std_block.path <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/covar_influence/"
uk.std_block.path <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/conservative_1e-5/covar_influence/"

#### First try: K=41
#Read in finngen
load(paste0(fg.std_block.path, "finngen_benchmark_svd_k41finngen_benchmark_svd_k41_unscaled_final_dat.RData"))
finngen.svd.41 <- scaleBoth(ret)
#Read in UKBB
load(paste0(uk.std_block.path, "ukbb_benchmark_svd_k41ukbb_benchmark_svd_k41_unscaled_final_dat.RData"))
ukbb.svd.41 <- scaleBoth(ret)
#Note- because its SVD, the columns should all be unit norm. 
#There are small numerical differences, so we just did normalized for goo dmeasure.
#stopifnot(max(unitNorms(ukbb.svd.41$V) - ukbb.svd.41$V) < 1e-12)
k41.unscaled <- alignAndCompare(finngen.svd.41,ukbb.svd.41)
cat("Correlation between V, K=41, unscaled setting:",round(k41.unscaled$adj.v, digits=3), "\n")
cat("Correlation between U, K=41, unscaled setting:",round(k41.unscaled$adj.u, digits=3), "\n")


#RRMSE:
rrmse.41.k <- getRRMSE(finngen.svd.41,ukbb.svd.41)
cat("RRMSE between X, K=41, unscaled setting:", round(rrmse.41.k, digits=3), "\n")

#### Try the scaled version- only slightly better
load(paste0(fg.std_block.path, "finngen_benchmark_svd_k41finngen_benchmark_svd_k41_scaled_final_dat.RData"))
finngen.svd.41.scaled <-scaleBoth(ret)
load(paste0(uk.std_block.path, "ukbb_benchmark_svd_k41ukbb_benchmark_svd_k41_scaled_final_dat.RData"))
ukbb.svd.41.scaled <- scaleBoth(ret)
res.scaled <- alignAndCompare(finngen.svd.41.scaled,ukbb.svd.41.scaled)
rrmse.41.res.scaled <- getRRMSE(finngen.svd.41.scaled,ukbb.svd.41.scaled)
cat("Correlation between V, K=41, scaled setting:",round(res.scaled$adj.v, digits=3), "\n")
cat("Correlation between U, K=41, scaled setting:",round(res.scaled$adj.u, digits=3), "\n")
cat("RRMSE between X, K=41, unscaled setting:", round(rrmse.41.res.scaled, digits=3), "\n")
#Only slightly better

#Try the K=19 vs K=17 run for a better comparison:
load(paste0(fg.std_block.path, "finngen_benchmark_svd_k19finngen_benchmark_svd_k19_unscaled_final_dat.RData"))
finngen.svd.19 <- scaleBoth(ret)
load(paste0(uk.std_block.path, "ukbb_benchmark_svd_k17ukbb_benchmark_svd_k17_unscaled_final_dat.RData"))
ukbb.svd.17 <- scaleBoth(ret)
res.19_v_17 <- alignAndCompare(finngen.svd.19,ukbb.svd.17)
rrmse.19_v_17 <- getRRMSE(finngen.svd.19,ukbb.svd.17)
cat("Correlation between V, K=19/17, UNscaled setting:",round(res.19_v_17$adj.v, digits=3), "\n")
cat("Correlation between U, K=19/17, UNscaled setting:",round(res.19_v_17$adj.u, digits=3), "\n")
cat("RRMSE between X, K=19/17, unscaled setting:", round(rrmse.19_v_17, digits=3), "\n")

#Try the K=19 vs K=17 scaled
load(paste0(fg.std_block.path, "finngen_benchmark_svd_k19finngen_benchmark_svd_k19_scaled_final_dat.RData"))
finngen.svd.19.scaled <- scaleBoth(ret)
load(paste0(uk.std_block.path, "ukbb_benchmark_svd_k17ukbb_benchmark_svd_k17_scaled_final_dat.RData"))
ukbb.svd.17.scaled <- scaleBoth(ret)
res.19_v_17.scaled <- alignAndCompare(finngen.svd.19.scaled,ukbb.svd.17.scaled)
rrmse.19_v_17.scaled <- getRRMSE(finngen.svd.19.scaled,ukbb.svd.17.scaled)
cat("Correlation between V, K=19/17, scaled setting:",round(res.19_v_17.scaled$adj.v, digits=3), "\n")
cat("Correlation between U, K=19/17, scaled setting:",round(res.19_v_17.scaled$adj.u, digits=3), "\n")
cat("RRMSE between X, K=19/17, scaled setting:", round(rrmse.19_v_17.scaled, digits=3), "\n")

########### This is the one with the lowest RRSME, so we go with this one.
## Most generous case- both are K=17
load(paste0(fg.std_block.path, "finngen_benchmark_svd_k17finngen_benchmark_svd_k17_unscaled_final_dat.RData"))
finngen.svd.17 <- scaleBoth(ret)
res.17_v_17 <- alignAndCompare(finngen.svd.17,ukbb.svd.17)
rrmse.17_v_17 <- getRRMSE(finngen.svd.17,ukbb.svd.17)
cat("Correlation between V, K=17/17, UNscaled setting:",round(res.17_v_17$adj.v, digits=3), "\n")
cat("Correlation between U, K=17/17, UNscaled setting:",round(res.17_v_17$adj.u, digits=3), "\n")
cat("RRMSE between X, K=17/17, UNscaled setting:", round(rrmse.17_v_17, digits=3), "\n")

#Try the K=17 vs K=17 scaled
load(paste0(fg.std_block.path, "finngen_benchmark_svd_k17finngen_benchmark_svd_k17_scaled_final_dat.RData"))
finngen.svd.17.scaled  <- scaleBoth(ret)
res.17_v_17.scaled <- alignAndCompare(finngen.svd.17.scaled,ukbb.svd.17.scaled)
rrmse.17_v_17.scaled <- getRRMSE(finngen.svd.17.scaled,ukbb.svd.17.scaled)
cat("Correlation between V, K=17/17, scaled setting:",round(res.17_v_17.scaled$adj.v, digits=3), "\n")
cat("Correlation between U, K=17/17, scaled setting:",round(res.17_v_17.scaled$adj.u, digits=3), "\n")
cat("RRMSE between X, K=17/17, scaled setting:", round(rrmse.17_v_17.scaled, digits=3), "\n")

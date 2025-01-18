## Script to simulate GWAS data from factors, loadings, and correlation structure.
#6 input files

#Returns the variance of the betas (overal variance)
#variance of true data of the total
#variance of the noise given the total.
#t is betas
#m is mean
#n is noise
varReport <- function(effects, true.signal, noise)
{
  t <- unlist(as.vector(effects))
  m <- unlist(as.vector(true.signal))
  n <- unlist(as.vector(noise))
  stopifnot(length(t) == length(n))
  #total variance, fraction from true signal, fraction from noise, var true /var noise, #var mu over sum var #var noise over sum var
  return(c(var(t), var(m)/var(t), var(n)/var(t), var(m)/var(n), var(m)/(var(n) +  var(m)),var(n)/(var(n) +  var(m)) ))
}

library("optparse")
library("data.table")
library("ggplot2")
library("magrittr")
option_list <- list(
  make_option(c("-i", "--input"), default="",
              help="Path to input yaml file"),
  make_option(c("-s", "--seed"), default=22,
              help="Set random seed, default is 22"),
  make_option(c("-o", "--output"), default = '', type = "character",
              help="Output dir")
)
#t = c("--input=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/udler_based_500/k4/noise2/udler2_realistic-high-1_r-75_noise2/udler2_realistic-high-1_r-75_noise2.yml")
t = c("--input=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/updated_variance_scaling/V6_U6_mafmixed_n100000.high_covar_1block_cont_scaling.yml",
      "--output=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/fast_runs//V6_U6_mafmixed_n100000.high_covar_1block_cont_scaling/fake")
t = c("--input=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/final_sims_june_2024/no_overlap/V101_U101_MAF-mix_eur_N-10000_RHO-none_No-none.yml",
      "--output=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/no_overlap/V101_U101_MAF-mix_eur_N-10000_RHO-none_No-none/sim1")
t = c("--input=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration//simulating_factors/custom_easy/yaml_files/final_sims_june_2024/special_2b_overlap//V105_U101_MAF-mix_eur_N-mixed_RHO-2b_mid_mixed_p_No-2b_high_no.yml",
      "--output=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/special_2b_overlap/V105_U101_MAF-mix_eur_N-mixed_RHO-2b_mid_mixed_p_No-2b_high_no//sim1")

#args <- parse_args(OptionParser(option_list=option_list), args = t)
args <- parse_args(OptionParser(option_list=option_list))
#set.seed(args$seed)
######################
cat("-----------------Generating simulated GWAS summary stats-------------------\n")
cat(" YAML file: ",args$input,"\n" )
cat(" Output destination", args$output,"\n" )

#read in relevant things
#Maybe a better way would be to specify a parameter file, that has all of this information in it.
#This way requires duplicates of each file in every simulation directory, which isn't what we want, is it.
#Specify a setup file.
yml <- read.table(args$input,header = FALSE,sep = ",") %>% set_colnames(c("n", "p"))
if("sim_ref" %in% yml$n)
{
  message("Simulation path data already provided; Will skip simulation generation")
  quit()
}




n_o <- as.matrix(fread(unlist(yml[which(yml$n == "samp_overlap"),2])))

#new.seed=args$seed + as.numeric(n_o[1,1]) #seed differs depending on the sample overlap matrix that is passed in and the iteration.
#SEed wasn't being unique
#Needed to update so more likely to be unique- based on the seed in the simulation
new.seed = as.numeric(paste0(as.character(args$seed), yml[which(yml$n == "seed_init"),2]))

#message("New sim version: seed set by setting + n: ", new.seed)
message("Constructing new simulation with seed ", new.seed)
set.seed(new.seed)

rho <- as.matrix(fread(unlist(yml[which(yml$n == "pheno_corr"),2]))) #Correlation between phenotypes
f <- as.matrix(fread(unlist(yml[which(yml$n == "factors"),2])))
l <-  as.matrix(fread(unlist(yml[which(yml$n == "loadings"),2]))) #gut says this should be more dense, huh.
noise.scaler = 1
if("noise_scaler" %in% yml$n)
{
	noise.scaler = as.numeric(yml[which(yml$n == "noise_scaler"),2])
	if(noise.scaler != 1)
	{message("WARNING: noise scaling is being modified......")}
}
#message(paste0("Noise scaler is ", noise.scaler))
maf.mat <- as.matrix(fread(unlist(yml[which(yml$n == "maf"),2])))
M <- nrow(n_o)
N <- nrow(l)

if("herit_scaler" %in% yml$n)
{
  #Has to correspond to the number of features
  ntraits = nrow(f)
  #message("Scaling X according to empirical trait heritabilities")
  herit.settings <- unlist(yml[which(yml$n == "herit_scaler"),2])

  #Effect size distribution heritabilities come from 2018 text:
  #"Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits"
  #Supplementary Table 4
  dis.traits <- c("Alzheimers","Asthma","Coronary artery disease", "Type 2 diabetes","Crohns disease", "Inflammatory bowel disease",
                      "Ulcerative colitis","Rheumatoid arthritis")
  dis.herit <- c(173.7,400.6,134.6,198.4,380.5,214.3,301.3,174.6)*10^-5

  continuous_traits <- c("Age at menarche","BMI","Height","Hip circumference","Waist circumference","Waist-to-hip ratio","HDL cholesterol",
                         "LDL cholesterol","Total cholesterol","Triglycerides","Child birth weight","Childhood obesity","Neuroticism")

  cont.herit <- c(12.4,18.0,14.6,8.0,8.1,6.5,15.2,21.8,22.7,23.1,9.1,95.6,0.9)*10^-5

  #quick.test

  #choose the one based on the herit.settings
  herit.factors <- switch(herit.settings,
         "disease" = sample(dis.herit,size=ntraits,replace = TRUE),
         "mixed" = c(sample(dis.herit,size=floor(ntraits/2),replace = TRUE),sample(cont.herit,size=floor(ntraits/2),replace = TRUE)),
         "continuous" = sample(cont.herit,size=ntraits,replace = TRUE),
         "ukbb_real" = rep(variance.scale,5)[1:ntraits])
  scaling.d <- diag(herit.factors)
  #need to multiply each column of mu by  sqrt(scaling.d[i]/var(i))

  #TODO- implement how this D gets applied!


}


#Set up the necessary data for SE estimation
n.mat <- matrix(do.call("rbind", lapply(1:N, function(x) diag(n_o))), nrow = N, ncol = M)
var.maf <- 2*(maf.mat)*(1-maf.mat)
S <- 1/sqrt(var.maf *n.mat) #This is the standard error, assuming the variance of our phenotypes is scaled to 1 (reasonable)
#maf.mat <- matrix(rep(maf$maf[1:500],10), nrow= 500, ncol = 10)

#Quick sanity check: do our estimated S correspond with the actual S?
if("se" %in% yml$n)
{
  se <- as.matrix(fread(unlist(yml[which(yml$n == "se"),2])))
  s.p <- data.frame("est" = S[,1], "t"=se[,1])
  ggplot(s.p, aes(x = est, y = t)) + geom_point() + theme_classic(15) + ylab("Actual SE")+ xlab("Sim")
  suppressMessages(ggsave(filename = paste0(args$output, ".se_sim_check.png")))
}


#First, calculate the correlation effect matrix
#This should be symmetric...
#if its not, that's an error on the part of
#C is m x m
C <- matrix(NA, nrow = M, ncol = M)
for(i in 1:nrow(n_o))
{
  for(j in 1:ncol(n_o))
  {
    C[i,j] <- (n_o[i,j]/(sqrt(n_o[i,i]) * sqrt(n_o[j,j]))) * rho[i,j]
  }
}
stopifnot(isSymmetric(C))
library(matrixcalc)
if(!is.positive.definite(C))
{
  message("C is not PD, adjusting now...")
  library(corpcor)
  C.updated <- make.positive.definite(C, tol = 1e-5)
  #percent change
  prop.change <- norm(C- C.updated, type = "F")/norm(C, type = "F")
  C <- C.updated
}
#if(!all(diag(C) == 1))
if(isFALSE(all.equal(diag(C),1)))
{
  message("Warning- not all diagonal elements == 1. This simulation is likely to give some weird results.")
  print(diag(C))
  print(which(diag(C)!=1))
  print(1-diag(C))
  print(paste0(args$output, ".c_matrix.txt"))
}
#NOTE: could also implement in terms of the matrix normal, a single line. That would work too, but I
#think (?) would be the same

#Create an image with the correlation structure:
source("/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/gwas_decomp_ldsc/src/plot_functions.R")
rownames(C) <- paste0("T", 1:nrow(C)); colnames(C) <- paste0("T", 1:ncol(C))

o <- plotCorrelationHeatmap(C, typin = "None", title = "Generated matrix of spurious covariance structure.", show.nums = TRUE)
suppressMessages(ggsave(plot = o, filename = paste0(args$output, ".spurious_covariance.png")))
write.table(x = C, file = paste0(args$output, ".c_matrix.txt"), quote = FALSE, row.names = FALSE)

#Next, generate all the effect sizes and "true" sigmas.
suppressMessages(library(MASS))
betas <- matrix(NA, nrow = N, ncol = M)
out.ses <- matrix(NA, nrow = N, ncol = M)
var.dat <- matrix(NA, nrow = N, ncol = 6)

#reset the seed- not sure if the above matrix functions require random seed, but just in case:
global.noise <- c()
global.signal <- c()
set.seed(args$seed)
all.noise <- NULL
mu.tot <- l %*% t(f)
#OLD VERSION- not scaling heritabilities here, so this isn't saved.
for(s in 1:N) #N is number of SNPs
{
  mu <- f %*% l[s,]
  #Variance on beta is standard error squared.
  se.i <- (diag(S[s,]))
  sigma <- se.i %*% C %*% se.i

  noise <- (noise.scaler * mvtnorm::rmvnorm(1, sigma = sigma))
  all.noise <- rbind(all.noise, noise)
  global.noise <- c(global.noise, unlist(as.vector(noise)))
  global.signal <- c(global.signal, unlist(as.vector(mu)))
  #When compared to a rnorm version when indpendent samples, this is the same. so that's fine.
  betas[s,] <- mu + t(noise)

  #Write this out somehwere useful.
  var.dat[s,] <- varReport(betas[s,], mu, noise)

  out.ses[s, ] <- diag(sigma) #variance of individual betas; context is now missing. Seems a bit.. odd.
  #With this mind, we could just be using the sigma hats directly, don't need to do this little dance.
}

#Updated version, where we add all the noise at the end after scaling for heritabilities.
betas.alt <- mu.tot + all.noise
stopifnot(max(betas.alt - betas) < 1e-10)
#Now scale mu to match our distributional assumptions
if("herit_scaler" %in% yml$n)
{
  #May updates- only scaling the non-zero SNPs to match the desired distribution:
  #get the scaling factors
  mu.tot.var <- apply(mu.tot, 2, function(x) var(x[x!=0])) #scaling the variance of causal SNPs only
  mu.scaled <- mu.tot

  #Cleaner way to code this:
  A <- diag(sqrt(herit.factors/mu.tot.var))
  mu.scaled <- l %*% t(f) %*% A #entries which are 0 won't be affected, so it works out just fine.
  #stopifnot(all(mu.scaled == mu.scaled.alt))
  #for(i in 1:length(mu.tot.var)) {mu.scaled[(mu.scaled[,i] !=0),i] <- mu.scaled[(mu.scaled[,i] !=0),i] * sqrt(herit.factors[i]/mu.tot.var[i])} #scaling just the variants that are CAUSAL
  #sanity check
  new.vars <- apply(mu.scaled, 2, function(x) var(x[x != 0]))
  stopifnot(max(new.vars - herit.factors) < 1e-16)
  #Alternative check- post multiply by A

  betas = mu.scaled + all.noise
  #Some plots if doing this manually

  pb <- cor(betas); rownames(pb) = paste0("T", 1:10); colnames(pb) = paste0("T", 1:10)
  suppressMessages(ggsave(plot = plotCorrelationHeatmap(pb,typin = "None",title = "Correlation structure of scaled beta hats"),
                          filename = paste0(args$output, ".effect_size_estimate_correlation.png")))

  mu.tot <- mu.scaled
}


###Write out simulations.
library(magrittr)
out.beta <- data.frame("SNP" = paste0("rs", 1:N), round(betas, digits = 7)) %>% set_colnames(c("SNP", paste0("T", 1:M)))
out.beta_true <- data.frame("SNP" = paste0("rs", 1:N), round(mu.tot, digits = 7)) %>% set_colnames(c("SNP", paste0("T", 1:M)))
stopifnot(max(abs((betas-mu.scaled) - all.noise)) < 1e-15)
#TODOIPDATE
out.se <- data.frame("SNP" = paste0("rs", 1:N), round(out.ses, digits = 7)) %>% set_colnames(c("SNP", paste0("T", 1:M)))
out.var <- data.frame("SNP" = paste0("rs", 1:N), round(var.dat, digits = 7)) %>%
  set_colnames(c("SNP", "var_beta", "var_mu:var_beta","var_noise:var_beta", "var_mu:var_noise", "var_m:var_m+var_n", "var_n:var_m+var_n"))

##Added in to run factorGo- write out Z and N:
out.z <- data.frame("SNP" = paste0("rs", 1:N), round(betas/out.ses, digits = 7)) %>% set_colnames(c("SNP", paste0("T", 1:M)))
out.n <- data.frame("N"=colMeans(n.mat))
#total variance, fraction from true signal, fraction from noise, var true /var noise, #var mu over sum var #var noise over sum vara
#Write out output:
write.table(x = out.beta, file = paste0(args$output, ".effect_sizes.txt"), quote = FALSE, row.names = FALSE)
write.table(x = out.se, file = paste0(args$output, ".std_error.txt"), quote = FALSE, row.names = FALSE)
write.table(x=var.dat, file = paste0(args$output, ".variance_data.txt"), quote = FALSE, row.names = FALSE)
write.table(x = C, file = paste0(args$output, ".c_matrix.txt"), quote = FALSE, row.names = FALSE)
#write out empirical covar
write.table(x = cor(all.noise), file = paste0(args$output, ".empirical_covar.txt"), quote = FALSE, row.names = FALSE)
global.var.report <- varReport(betas, global.signal, global.noise)

#Write out for factorGo
write.table(x=out.z, file = paste0(args$output, ".z.txt"), quote = FALSE, row.names = FALSE,sep = "\t")
write.table(x=out.n, file = paste0(args$output, ".N.txt"), quote = FALSE, row.names = FALSE)

#Write global data, if its not there.....
global_beta = paste0(gsub(args$output, pattern = "sim\\d+", replacement = ""),"noise-free_effect_sizes.txt")
if(!file.exists(global_beta))
{
  #message("writing global beta without noise file.")
  write.table(x=out.beta_true, file = global_beta, quote = FALSE, row.names = FALSE)

}

sink(paste0(args$output, ".variance_notes.txt"))
cat(paste0("Globally, ", round(mean(global.var.report[5])*100, digits =2), "% of sample variation from true value."))
cat(paste0("Globally, ", round(mean(global.var.report[6])*100,digits = 2), "% of sample variation from noise."))
cat(paste0("Globally, signal to noise ratio is:", round(mean(global.var.report[4]),digits = 4)))
cat(paste0("On average, ", round(mean(var.dat[,5])*100, digits =2), "% of sample variation from true value.\n"))
cat(paste0("On average, ", round(mean(var.dat[,6])*100,digits = 2), "% of sample variation from noise."))
cat(paste0("On average, signal to noise ratio is:", round(mean(var.dat[,4]),digits = 4)))
sink()

message(paste0("Globally, ", round(mean(global.var.report[5])*100, digits =2), "% of sample variation from true value."))
message(paste0("Globally, ", round(mean(global.var.report[6])*100,digits = 2), "% of sample variation from noise."))
message(paste0("On average, ", round(mean(var.dat[,5])*100, digits =2), "% of sample variation from true value."))
message(paste0("On average, ", round(mean(var.dat[,6])*100,digits = 2), "% of sample variation from noise."))
#message(paste0("On average, signal to noise ratio is:", round(mean(var.dat[,4]),digits = 4)))
#need a good way to look at SEs....

#The 2nd PC correlates strongly with sample size
#the 1st PC correlates strongly with MAF?
cat("----------------------------------------------------------\n")


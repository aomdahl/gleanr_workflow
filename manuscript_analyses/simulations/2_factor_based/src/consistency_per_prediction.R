suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("magrittr"))
suppressMessages(library("dplyr"))
suppressMessages(library(stringr))
suppressMessages(library("combinat"))
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/matrix_similarity_scoring/evaluateSimR2.R")
library("gleanr")
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/matrix_similarity_scoring/matrix_comparison_utils.R")

#updated this to reflect cleaner code in https://stackoverflow.com/questions/68773192/how-to-calculate-variance-between-observations-when-observations-are-in-matrix
varAllEntries <- function(aligned_matrices, nrow, ncol, col_tag="T")
{
  #full_array <- array(unlist(aligned_matrices), c(nrow, ncol, length(aligned_matrices)))
  full_array <- do.call(abind::abind, c(aligned_matrices, along = 3))
  var.by.row.by.col <- apply(full_array, 1:2, var) #go along dimensions 1 and 2 to get the vector
  #This is equivilant to the following, but we opt for the shorter thing
  #var.by.row.by.col <- do.call("rbind", lapply(1:nrow(full_array), function(i) {
  #  sapply(1:ncol(full_array), function(j) var(full_array[i,j,]))
  #}))

  colnames(var.by.row.by.col) <- paste0("F", 1:ncol)
  rownames(var.by.row.by.col) <- paste0(col_tag, 1:nrow)
  var.by.row.by.col
}

pairwiseCorByMatrix <- function(aligned_matrices)
{
  full.len <- expand.grid(1:10, 1:10) %>% filter(Var1 != Var2)
  overall.r2 <- c()
  for(i in 1:nrow(full.len))
  {
    overall.r2 <- c(overall.r2, evaluateSingleMatrixConstruction(aligned_matrices[[full.len[i,1]]], aligned_matrices[[full.len[i,2]]]))
  }
  list("avg" = mean(overall.r2), "var"=var(overall.r2))
    
}


varByMatrix <- function(aligned_matrices)
{
  #First, get the "average" matrix
  full_array <- do.call(abind::abind, c(aligned_matrices, along = 3))
  mean.mat <- apply(full_array, 1:2, mean) #go along dimensions 1 and 2 to get the vector
  sum(sapply(aligned_matrices, function(x) norm(x-mean.mat, "F")^2))/length(aligned_matrices)
  #Updated to say squared, as it should!
  #sum(sapply(aligned_matrices, function(x) norm(x-mean.mat, "F")))/length(aligned_matrices)
}

varByFactor <- function(aligned_matrices)
{
  #First, get the "average" matrix
  full_array <- do.call(abind::abind, c(aligned_matrices, along = 3))
  mean.mat <- apply(full_array, 1:2, mean) #go along dimensions 1 and 2 to get the vector

  sapply(1:ncol(mean.mat), function(i)
    {
      sum(sapply(aligned_matrices, function(x) norm(matrix(x[,i]-mean.mat[,i]), "F")))/length(aligned_matrices)
    })
}


suppressMessages(library(NatParksPalettes))
option_list <- list(
  make_option(c("-o", "--output"), default = '', type = "character",
              help="Output dir + handle"),
  make_option(c("-s", "--sim_path"), default = '', type = "character",
              help="File containing all of the paths of simulations you wish to go into this analysis"),
  make_option(c("-y", "--yaml"), default = '', type = "character",
              help="Path to the settings file."),
  make_option(c("-n", "--scale_data"), default = FALSE, type = "logical",
              help="Scale true and predicted V/U to have unit norm.", action = "store_true"),
  make_option(c("--sim_raw"), default = "", type = "character",
              help="Path to the actual simulation files, not the tests")

)

#debug case
finalpath="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/1b_overlap/"
t = c(paste0("--output=", finalpath, "V101_U102_MAF-mix_eur_N-50000_RHO-1b_high_mixed_p_No-1b_high_no/factorization_results/summary"),
      "--yaml=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/final_sims_june_2024/1b_overlap/V101_U102_MAF-mix_eur_N-50000_RHO-1b_high_mixed_p_No-1b_high_no.yml",
      "--scale_data", paste0("--sim_path=", finalpath, "V101_U102_MAF-mix_eur_N-50000_RHO-1b_high_mixed_p_No-1b_high_no/factorization_results/"))

test_path="/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/final_sims_june_2024/special_mid-density_no_overlap/"
t = c(paste0("--output=", test_path, "V111_U103_MAF-mix_eur_N-5000_RHO-none_No-none/factorization_results/summary"),
      "--yaml=/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/yaml_files/final_sims_june_2024/special_mid-density_no_overlap/V111_U103_MAF-mix_eur_N-5000_RHO-none_No-none.yml",
      "--scale_data", paste0("--sim_path=", test_path, "V111_U103_MAF-mix_eur_N-5000_RHO-none_No-none//factorization_results/"))


#args <- parse_args(OptionParser(option_list=option_list), args = t)
######

args <- parse_args(OptionParser(option_list=option_list))#, args = t)

if(args$sim_raw == "")
{
  message("setting sim raw directory as one up from sim_path")
  args$sim_raw = gsub(pattern="factorization_results\\/",x=args$sim_path, replacement="")
}
yml <- read.table(args$yaml, sep = ",")
methods.run <- ((yml %>% filter(V1 == "test_methods"))$V2 %>% strsplit(.,":" ))[[1]]
true.loadings <- as.matrix(fread((yml %>% filter(V1 == "loadings"))$V2))
true.factors <-  as.matrix(fread((yml %>% filter(V1 == "factors"))$V2))

#Scale data, should be always true b/c
if(args$scale_data)
{
  true.loadings <- unitNorms(true.loadings)
  true.factors <- unitNorms(true.factors)
}
niter <- as.numeric((yml %>% filter(V1 == "iter"))$V2)
s <- args$sim_path
true.k <- ncol(true.factors)
#This will evaluate the simulations you desire. The assumption is that each individual simulation has its own directory, with a
#sim_loadings and sim_factors file, and then files corresponding to the predictved outcome for several different methods.
sim.performance <- NULL
f_i =1
nmethods = length(methods.run)
ref.list <- list()
pred.list <- list()

#Output data storage:
factorVar <- NULL;matrixVar<- NULL;perCellVar<- NULL;perMatAvg <- NULL;

n_features=12
for(m in methods.run)
{
  message("Assessing simulations from ", m)
  r_performance <- matrix(NA, nrow = niter, ncol = n_features)
  all_iter_aligned <- list()
  for(i in 1:niter){
    lookup_id=paste0(m,"_", i)
    pred.list[[lookup_id]] <- list()
    ref.list[[lookup_id]] <- list("U"=true.loadings, "V"=true.factors)
    if(!file.exists(paste0(s, "/sim",i, ".", m, ".loadings.txt")) & i > 1)
    {
      stop("WARNING: missing file: ", paste0(s, "/sim",i, ".", m, ".loadings.txt"))
      next
    }
    if(!file.exists(paste0(s, "/sim",i, ".", m, ".loadings.txt")) & i == 1)
    {
      stop("WARNING: missing file: ",paste0(s, "/sim",i, ".", m, ".loadings.txt"))
      break;
    }
    pred.list[[lookup_id]][["U"]] <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".loadings.txt")))
    pred.list[[lookup_id]][["V"]] <- as.matrix(fread(paste0(s, "/sim",i, ".", m, ".factors.txt")))
    if(args$scale_data)
    {
      # message("Scaling both true and loaded data for convenient comparison")
      if(!all(pred.list[[lookup_id]][["U"]] == 0))
      {
        pred.list[[lookup_id]][["U"]] <-unitNorms(pred.list[[lookup_id]][["U"]])
        pred.list[[lookup_id]][["V"]] <- unitNorms(pred.list[[lookup_id]][["V"]])
      }else
      {
        message("Unable to scale data; all elements 0'd out.")
      }

    }
    all_iter_aligned[[i]] <- evaluate_error(true.loadings, true.factors, pred.list[[lookup_id]][["U"]], pred.list[[lookup_id]][["V"]])
    f_i = f_i+1
  }
  
  #JJune 21, added in comments to try and clean things up for faster running
  ### Calculate the consistency
  v.matrices <- lapply(all_iter_aligned, function(x) x$reordered_V)
  u.matrices <- lapply(all_iter_aligned, function(x) x$reordered_U)
  v.vars.all <- varAllEntries(v.matrices,nrow(true.factors), ncol(true.factors), col_tag="Trait")
  u.vars.all <- varAllEntries(u.matrices,nrow(true.loadings), ncol(true.loadings), col_tag="SNP")
  #r2.avg.v <- pairwiseCorByMatrix(v.matrices) 
  #r2.avg.u <- pairwiseCorByMatrix(u.matrices)
  ##Save this data in a cohesive way
  factorVar <- rbind(factorVar, data.frame("method"=m, "matrix"="V","factor_vars"= varByFactor(v.matrices)),
                     data.frame("method"=m, "matrix"="U","factor_vars"=varByFactor(u.matrices) ))

  matrixVar <- rbind(matrixVar, data.frame("method"=m, "matrix"="V","matrix_vars"= varByMatrix(v.matrices)),
                     data.frame("method"=m, "matrix"="U","matrix_vars"=varByMatrix(u.matrices) ))


  #perCellVar <- rbind(perCellVar,
  #                    data.frame("method"=m, "matrix"="V",reshape2::melt(v.vars.all)),
  #                    data.frame("method"=m, "matrix"="U",reshape2::melt(u.vars.all)))
  
  
  #perMatAvg <- rbind(perMatAvg,
  #                    data.frame("method"=m, "matrix"="V","mean"=r2.avg.v$avg, "var"=r2.avg.v$var),
  #                    data.frame("method"=m, "matrix"="U","mean"=r2.avg.u$avg, "var"=r2.avg.u$var))

}

#Clean up tables if needed
stopifnot(ncol(perCellVar) == 5)
#perCellVar %<>% set_colnames(c("method","matrix", "Trait", "Factor","entry_var"))
#Now plot it, if desired
#write.table(factorVar, file = paste0(args$output, ".per_factor_variance.tsv"), quote = FALSE, row.names = FALSE)
write.table(matrixVar, file = paste0(args$output, ".per_matrix_variance.tsv"), quote = FALSE, row.names = FALSE)
#write.table(perCellVar, file = paste0(args$output, ".per_entry_variance.tsv"), quote = FALSE, row.names = FALSE)
#write.table(perMatAvg, file = paste0(args$output, ".per_matrix_pairwise_r2.tsv"), quote = FALSE, row.names = FALSE)
##For interactive use
if(FALSE)
{
  ggplot(factorVar, aes(x=method, y = factor_vars, fill=method)) + geom_boxplot() +
    facet_wrap(~matrix) + theme_bw() +
    xlab("Factorization method")+ ylab("Per-factor variance across simulations")+
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))

  ggplot(perCellVar, aes(x=method, y = entry_var, fill=method)) + geom_boxplot() +
    facet_wrap(~matrix) + theme_bw() +
    xlab("Factorization method")+ ylab("Per-cell variance across simulations")+
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))

  ggplot(matrixVar, aes(x=method, y = matrix_vars, fill=method)) + geom_bar(stat="identity") +
    facet_wrap(~matrix) + theme_bw() +
    xlab("Factorization method")+ ylab("Per-cell variance across simulations")+
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))

}





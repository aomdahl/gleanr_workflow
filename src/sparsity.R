
exactSparsity <- function(dat){
  stats <- colSums(dat == 0)/nrow(dat)
  #summary(stats)
  #plot(stats)  
  return(stats)
}

approxSparsity <- function(dat, thresh = 1e-4){
  stats <- colSums(abs(dat) < thresh)/nrow(dat)
  return(stats)
}

exactAvgSparsity <- function(dat){
  sum(colSums(dat == 0))/nrow(dat) /ncol(dat)
}

approxSparsity <- function(dat, thresh = 1e-4){
  stats <- colSums(abs(dat) < thresh)/nrow(dat)
  return(stats)
}

approxAvgSparsity <- function(dat, thres){
  sum(colSums(abs(dat) < thresh))/nrow(dat) /ncol(dat)
}
#quick script to choose the traits we want
pacman::p_load(data.table, tidyr, dplyr,  ggplot2, stringr, argparse)
parser <- ArgumentParser$new()
parser$add_description("Script to quickly assess the sparsity of a matrix by column")
parser$add_argument("--data", type = 'character', help = "Path to data file")
parser$add_argument('-hist',type='character', help='Do you want to make a figure with the data? give the path', default = "")
parser$add_argument('--figure',type='character', help='Do you want to make a figure with the data? give the path', default = "")
parser$add_argument('--approxthresh',type='numeric',help='Specify the sparsity threshold', default = 1e-4)
parser$add_argument('--approx',type='logical',action='store_true', help='Give approximate sparsity', default = FALSE)
parser$add_argument('--omit_first',type='logical',action='store_true', help='Omit the first column in calculation', default = FALSE)
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
#options(readr.num_columns = 0)

if(FALSE)
{
  args <- list()
  args$data <- "/Users/ashton/Documents/JHU/Research/LocalData/snp_network/Exploratory_notebooks/alternative_factors/sparsePCA_alpha0.001.txt"
  args$figure <- TRUE
}

dat <- fread(args$data)
if(args$omit_first)
{
  dat <- dat[,-1]
}
if(args$approx)
{
  #Approximate sparsity
  sparsity <- approxSparsity(dat, thresh = args$approxthresh)
  print("Approximate sparsity:")
  print(sparsity)
} else{
  #Exact sparsity 
  sparsity <- exactSparsity(dat)
  print("Exact sparsity:")
}

print(summary(sparsity))
if(args$figure != "")
{
    png(args$figure)
    plot(sparsity)
    dev.off()

}
if(args$hist != "")
{
    png(args$hist)
    hist(as.vector(as.matrix(dat)), main = "Distribution of values", breaks = 50, xlab = "Matrix entries")
    dev.off()

}

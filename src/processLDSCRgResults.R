pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr,stats, cowplot, magrittr, argparse)

renameAll <- function(m, n)
{
  colnames(m) <- n
  rownames(m) <- n
  m
}


parser <- ArgumentParser$new()
parser$add_description("Script to process the output of LDSC")
parser$add_argument("--input_path", type = 'character', help = "Path to tabular output")
parser$add_argument("--which_data", type = 'character', help = "Which data to plot and extract; default is all.")
parser$add_argument("--focus_trait", type = 'character', help = "specify the name of a focus trait, if you want", default = "")
parser$add_argument("--cohort_data", type = 'character', help = "specify the cohort of the traits, if you want", default = "")
parser$add_argument("--output", type = 'character', help = "Specify the path to the genomic correction coefficients. If none provided, none used", default = "")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/postMungeUtils.R")
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc//src/plot_functions.R")
if(TRUE)
{
  args <- list()
  args$output <- "/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/infertility_analysis/hm3_ldsc/tabular/TEST"
  args$input_path <- "/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/infertility_analysis/hm3_ldsc/tabular/"
  args$trait_names <- "default"
  args$which_data <- "all"
  args$focus_trait <- "female_infertility_ICD10"
  args$cohort_data <-"/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/infertility_analysis/trait_data/simple.second.cohort_data.txt"
    
  
}

#TODO: check that rows and columns are in the same order.
message("Evaluating genetic correlation")
#MAke the r2g and heatmap....

if(args$which_data == "intercept")
{
  intercepts <- ldscStdIntoMatrix(look_path = args$input_path,"" , filter_se = TRUE)
}else if(args$which_data == "rg")
{
  rg <- ldscGCOVIntoMatrix(look_path = args$input_path,"rg")
  if(args$trait_names == "default")
  {
    cn <- unlist(lapply(str_split(colnames(rg),pattern = "\\."), function(x) x[(length(x) - 2)]))
    rg <- renameAll(rg, cn)
  }
  plotCorrelationHeatmap(rg,typin ="heatmap.default", title = "Heatmap of genetic correlation (rg)")
  ggsave(filename = paste0(args$output, "rg_heatmap.png"),width = 10, height = 10)
  plotCorrelationHeatmap(rg,typin ="heatmap.default", title = "Heatmap of genetic correlation (rg)",p.vals = p)
  ggsave(filename = paste0(args$output, "rg_heatmap_with_p.png"),width = 15, height = 10)
  
  
} else if(args$which_data == "p")
{
  p <- ldscGCOVIntoMatrix(look_path = args$input_path,"p")
  if(args$trait_names == "default")
  {
    cn <- unlist(lapply(str_split(colnames(rg),pattern = "\\."), function(x) x[(length(x) - 2)]))
    p <- renameAll(p, cn)
  }
  
} else if(args$which_data == "h2_obs")
{
  h2obs <- ldscGCOVIntoMatrix(look_path = args$input_path,"h2_obs")
} else {
  mesasge("doing all by default")
  h2obs <- ldscGCOVIntoMatrix(look_path = args$input_path,"h2_obs")
  p <- ldscGCOVIntoMatrix(look_path = args$input_path,"p")
  rg <- ldscGCOVIntoMatrix(look_path = args$input_path,"rg")
  intercepts <- ldscStdIntoMatrix(look_path = args$input_path,"" , filter_se = TRUE)
  if(args$trait_names == "default")
  {
    cn <- unlist(lapply(str_split(colnames(rg),pattern = "\\."), function(x) x[(length(x) - 2)]))
    rg <- renameAll(rg, cn)
    p <- renameAll(p, cn)
  }
  
}





if(args$focus_trait != "")
{
  focus_cor <- data.frame("traits" = rownames(rg), "rg" = rg[,args$focus_trait], "p"= p[,args$focus_trait])
  ggplot(focus_cor, aes(x = reorder(traits, rg), y=rg, fill = -log10(p))) + geom_bar(stat = "identity") + theme_classic(15) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + xlab("GWAS traits") + ylab(expression(r[g])) + 
    ggtitle(paste0("Trait correlation with ", args$focus_trait)) + scale_fill_gradient(limits = c(0,max(1.3, max(-log10(focus_cor$p))))) + coord_flip()
  ggsave(filename = paste0(args$output, "correlation_", args$focus_trait, ".png"), width = 8, height = 10)
}

if(args$which_data == "ALL" | args$which_data == "h2_obs")
{
  
  ggplot(h2obs, aes(x = reorder(Trait, -h2_obs), y = h2_obs)) + geom_bar(stat ="identity") +
    geom_errorbar(aes(ymin=h2_obs-1.96*h2_obs_se, ymax=h2_obs+1.96*h2_obs_se), width=.2, position=position_dodge(.9)) + theme_classic() + coord_flip() + 
    xlab("Trait") + ylab("h2 (observed scale)")
  
  
  p<- ggplot(df2, aes(x=dose, y=len, fill=supp)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                  position=position_dodge(.9)) 
  
}

message("Now looking at cohort overlaps...")
int <- ldscGCOVIntoMatrix(look_path = args$input_path,"gcov_int",filter_se = TRUE)
int <- renameAll(int, cn)
if(args$cohort_data != "")
{
  cohorts <- fread(args$cohort_data) %>% set_colnames(c("trait_name", "Cohort")) %>% 
    mutate("trait_ordered" = factor(trait_name, levels = cn)) %>% drop_na() %>% arrange(trait_ordered)
}
plotCorrelationHeatmap(int,typin ="heatmap.default", title = "Heatmap of LDSC Rg Intercept\n(cohort overlap effects)",colors = cohorts$Cohort)
ggsave(filename = paste0(args$output, "gcov_int.png"),width = 15, height = 10)

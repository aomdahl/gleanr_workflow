pacman::p_load(data.table, tidyr, dplyr, ggplot2, stringr,stats, cowplot, magrittr, optparse)

renameAll <- function(m, n)
{
  colnames(m) <- n
  rownames(m) <- n
  m
}


option_list <- list(
make_option(c("--input_path"), type = 'character', help = "Path to tabular output"),
make_option(c("--which_data"), type = 'character', help = "Which data to plot and extract; default is all.", default = "ALL"),
make_option(c("--focus_trait"), type = 'character', help = "specify the name of a focus trait, if you want", default = ""),
make_option(c("--cohort_data"), type = 'character', help = "specify the cohort of the traits, if you want", default = ""),
make_option(c("--output"), type = 'character', help = "Get the specific data you're interested in...", default = ""),
make_option(c("--filter_se"), type = 'logical', action = "store_true", help = "Specify this if you want to filter your outputs by the standard error.", default = FALSE)
)
args <- parse_args(OptionParser(option_list=option_list))
args$trait_names <- "default" #haven't implemented this yet...
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/postMungeUtils.R")
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc//src/plot_functions.R")
if(FALSE)
{
  args <- list()
  args$output <- "/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/infertility_analysis/hm3_ldsc/tabular/TEST"
  args$input_path <- "/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/infertility_analysis/hm3_ldsc/tabular/"
  args$trait_names <- "default"
  args$which_data <- "gcov_int"
  args$focus_trait <- "female_infertility_ICD10"
  args$cohort_data <-"/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/infertility_analysis/trait_data/simple.second.cohort_data.txt"
  args$filter_se <- FALSE
  args$trait_names <- "default" #haven't implemented this yet..
  
}

if(FALSE)
{
  args <- list()
  args$output <- "/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/custom_l1_factorization/results/udler_original/ldsc_summary"
  args$input_path <- "/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/udler_td2/ldsr_all/ldsr_results/"
  args$trait_names <- "default"
  args$which_data <- "gcov_int"
  args$cohort_data <-"/home/aomdahl1/scratch16-abattle4/ashton/snp_networks/scratch/infertility_analysis/trait_data/simple.second.cohort_data.txt"
  args$filter_se <- FALSE
  args$trait_names <- "default" #haven't implemented this yet..
}

#TODO: check that rows and columns are in the same order.

#MAke the r2g and heatmap....

if(args$which_data[1] == "intercept" | args$which_data == "ALL")
{
  message("looking for intercept")
  intercepts <- ldscStdIntoMatrix(look_path = args$input_path,"" , filter_se = args$filter_se)
  write.table(intercepts, file = paste0(args$output, "genomic_correction_intercept.tsv"), quote = FALSE, row.names = FALSE)
}
if(args$which_data == "rg"| args$which_data == "ALL")
{
  rg <- ldscGCOVIntoMatrix(look_path = args$input_path,"rg",filter_se = args$filter_se)
  if(args$trait_names == "default")
  {
    cn <- unlist(lapply(str_split(colnames(rg),pattern = "\\."), function(x) x[(length(x) - 2)]))
    rg <- renameAll(rg, cn)
  }
  plotCorrelationHeatmap(rg,typin ="heatmap.default", title = "Heatmap of genetic correlation (rg)")
  ggsave(filename = paste0(args$output, "rg_heatmap.png"),width = 10, height = 10)
  
}
if(args$which_data == "p"| args$which_data == "ALL")
{
  p <- ldscGCOVIntoMatrix(look_path = args$input_path,"p")
  if(args$trait_names == "default")
  {
    cn <- unlist(lapply(str_split(colnames(rg),pattern = "\\."), function(x) x[(length(x) - 2)]))
    p <- renameAll(p, cn)
  }
  if(args$which_data == "ALL") {
    diag(p) <- 0
  plotCorrelationHeatmap(rg,typin ="heatmap.default", title = "Heatmap of genetic correlation (rg)",p.vals = p) + ggplot2::scale_alpha_continuous(guide = "none")
  #se <- ldscGCOVIntoMatrix(look_path = args$input_path,which.col = "rg_se")
    #corrplot::corrplot(cormat_mod, method = "square", p.mat = p.vals[order_dat$order_x, order_dat$order_y]* 0.0001)
    # ldscTableToMatrix <- function(joined_main, col_pick, diag_scores = 1, null_scores = 0)
    #corrplot::corrplot(cormat_mod,lowCI.mat = cormat_mod -se, uppCI.mat = cormat_mod + se )
    
    
  ggsave(filename = paste0(args$output, "rg_heatmap_with_p.png"),width = 15, height = 10)
  }
}
if(args$which_data == "h2_obs"| args$which_data == "ALL")
{
  h2obs <- ldscGCOVIntoMatrix(look_path = args$input_path,"h2_obs")
} else {
  message("doing all by default")
  h2obs <- ldscGCOVIntoMatrix(look_path = args$input_path,"h2_obs")
  p <- ldscGCOVIntoMatrix(look_path = args$input_path,"p")
  rg <- ldscGCOVIntoMatrix(look_path = args$input_path,"rg", filter_se = args$filter_se,filter_fdr = 0.05 )
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
  ggsave(filename = paste0(args$output, "_heritability_obs.png"), width = 8, height = 10)
  
  
 # p<- ggplot(df2, aes(x=dose, y=len, fill=supp)) + 
 #   geom_bar(stat="identity", color="black", 
 #            position=position_dodge()) +
 #   geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
 #                 position=position_dodge(.9)) 
  
}

if(args$which_data == "gcov_int"| args$which_data == "ALL")
{
  message("Now looking at cohort overlaps...")
  int <- ldscGCOVIntoMatrix(look_path = args$input_path,"gcov_int",filter_se = args$filter_se, filter_fdr = 0.05)
  if(args$trait_names == "default")
  {
    cn <- unlist(lapply(str_split(colnames(int),pattern = "\\."), function(x) x[(length(x) - 2)]))
    int <- renameAll(int, cn)
  }
print(cn) 
  if(args$cohort_data != "")
  {
    cohorts <- fread(args$cohort_data) %>% set_colnames(c("trait_name", "Cohort")) %>% 
      mutate("trait_ordered" = factor(trait_name, levels = cn)) %>% drop_na() %>% arrange(trait_ordered)
    #make the colors tha tmatch too...
    library(RColorBrewer)
    num <- length(unique(cohorts$Cohort))
    if(num <=12)
    {
      col <- brewer.pal(num, "Paired")
    }
    merge.frame <- data.frame("col" = col, "Cohort" = unique(cohorts$Cohort))
    cohorts <- left_join(cohorts, merge.frame, by = "Cohort")
  plotCorrelationHeatmap(int, title = "Heatmap of LDSC Rg Intercept\n(cohort overlap effects)",colors = cohorts$col)
  ggsave(filename = paste0(args$output, "gcov_int.png"),width = 15, height = 10)
  

  }else
  {
  print(int)
	  plotCorrelationHeatmap(int, title = "Heatmap of LDSC Rg Intercept\n(cohort overlap effects)")
  
  ggsave(filename = paste0(args$output, "gcov_int.png"),width = 15, height = 10)
  
  }
  #typin ="heatmap.default"
  
  
  #Save also the matrix itself.
  stopifnot(all(colnames(int) == rownames(int)))
  write.table(int, file= paste0(args$output, "gcov_int.tab.csv"), quote = FALSE,row.names = FALSE)
}

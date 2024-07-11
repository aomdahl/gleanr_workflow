####### Visualize LDSC Enrichments
####### Made by Ashton Omdahl, Battle Lab, 2020-2021
#######
#deal with labels, add in GTEX colors, etc.
labelsPrep <- function()
{
  labels <- fread("/data/abattle4/aomdahl1/reference_data/finucane_2018_supp_table7.ashton.csv")
  color_labels <- fread("/data/abattle4/aomdahl1/reference_data/finucane_2018_supp_table7.tissue_color_assigns.csv")
  tissue_names <- fread("/data/abattle4/aomdahl1/reference_data/finucane_2018_supp_table7.tissue_names_pretty.csv")
  #Kyped the following code from Rebecca Keener on 6/2/2021
  yes<-subset(labels, labels$entex=="Yes")
  yes$Name<-paste(yes$tissue, "_ENTEX__", yes$mark, sep="")
  no<-subset(labels, labels$entex=="No")
  no$Name<-paste(no$tissue, no$mark, sep="__")
  labels<-rbind(yes, no)
  labels <- left_join(labels, tissue_names, by = "tissue") %>% left_join(., color_labels, by = "new_category")
  labels
}



pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, optparse)
#"Script to visualize LDSC-cell specific results"
option_list <- list(
make_option("--extension", type = 'character', help = "File extension regex", default ="*.multi_tissue.cell_type_results.txt" ),
make_option("--output", type = 'character', help = "output file to write to"),
make_option("--input_file", type = 'character', help = "Input file- typically used for baseline plots"),
make_option("--input_dir", type = 'character', default = "/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc/results/ldsc/seed1_thresh0.7/", help = "where to look"),
make_option("--plot_type", type = 'character', default = "vertical", help = "What kind of plot would you like?"),
make_option("--tissue_db", type = 'character', default = "Multi_tissue_chromatin", help = "Source did you use for the tissue-specific annotations? Multi_tissue_gene_expr or chromatin_marks"),
make_option("--fdr", type = 'numeric', default = 0.05, help = "FDR threshold for trait enrichment"),
make_option("--orientation", type = 'character', default = "vertical", help= "Specify if you want a vertical or horizontal plot."),
make_option("--include_factors", type = "character", default = "ALL", help = "Specify which factors to include, comma separated, if only interested in some. Note that this does affect FDR calculation.")
)

test=c("--output=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/NONE_ldsc_enrichment_Multi_tissue_chromatin/",
       "--input_dir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/NONE_ldsc_enrichment_Multi_tissue_chromatin/",
       "--plot_type=data_tab")

args <- parse_args(OptionParser(option_list=option_list))

dir <- args$input_dir
ext <- args$extension

message("Code needs to be cleaned-up")
substitution_regex = "__"
substitution_regex_to = "-"
if(args$plot_type == "baseline")
{
    ldsc_output <- fread(args$input_file)
    ldsc_output$padj <- p.adjust(ldsc_output$Enrichment_p, method = "bonferroni")
    ldsc_output$logpadj <- -log(ldsc_output$padj)
    ldsc_output$names <- gsub("_0", "", ldsc_output$Category)
    p <- ggplot(ldsc_output %>% drop_na() %>% filter(padj != 1), aes(x = reorder(names, Enrichment), y = Enrichment, fill = logpadj)) + geom_bar(stat = "identity")  + theme_minimal(15) + ggtitle("Enrichment") + xlab("Baseline Annotation") + ylab("Enrichment score") + labs(fill = "Enrichment \n adjusted \n -log10(p)") + scale_fill_gradient(low = "gray97", high = "red")+ coord_flip()
    ggsave(args$output, plot = p, width = 7, height = 9)
    quit()
}

file.list <- list.files(dir, pattern = ext)
if(length(file.list) == 0)
{
    print(file.list)
    print(dir)
    print(ext)
    print("The path you provided didn't quite work. Please check it and try again")
    quit()
}

all <- lapply(file.list, function(x) fread(paste0(dir, x)) %>% arrange("Name") %>% 
                mutate("Source" = str_extract(str_split(x, pattern = "\\.")[[1]][1], pattern = "F\\d+")))  

#Where do the labels come from? Supplementary figure in the paper,
if(args$tissue_db == "Multi_tissue_gene_expr")
    {
  #Importantly- fix the names so they all match up  
  
todos <- plyr::rbind.fill(all) %>% mutate("-log10(P)" =-log10(Coefficient_P_value)) %>% 
  mutate("Name" = gsub(pattern = "_",replacement = " ",x = Name))
      #https://www.nature.com/articles/s41588-018-0081-4#Sec40
        #labels <- fread("/scratch16/abattle4/ashton/snp_networks/ldsc_tissue_type_reference.txt") %>% rename("Name" = Tissue)
        labels <- fread("/scratch16/abattle4/ashton/snp_networks/ldsc_tissue_type_reference.ashton.txt") %>% rename("Name" = tissue)
        #Optional name mod:
        labels$source <- ifelse(grepl("A\\d+\\.", labels$Name), "(Franke)", "(GTEx)")
        labels$Name_source <- paste0(labels$Name, labels$source)
        substitution_regex = "A\\d+\\.(\\d\\d\\d.)*"
        substitution_regex_to = ""
    } else {
      todos <- plyr::rbind.fill(all) %>% mutate("-log10(P)" =-log10(Coefficient_P_value))
      if(FALSE)
      {
        
        labels <- fread("/data/abattle4/aomdahl1/reference_data/finucane_2018_supp_table7.ashton.csv")
        color_labels <- fread("/data/abattle4/aomdahl1/reference_data/finucane_2018_supp_table7.tissue_color_assigns.csv")
        tissue_names <- fread("/data/abattle4/aomdahl1/reference_data/finucane_2018_supp_table7.tissue_names_pretty.csv")
        #Kyped the following code from Rebecca Keener on 6/2/2021
        yes<-subset(labels, labels$entex=="Yes")
        yes$Name<-paste(yes$tissue, "_ENTEX__", yes$mark, sep="")
        no<-subset(labels, labels$entex=="No")
        no$Name<-paste(no$tissue, no$mark, sep="__")
        labels<-rbind(yes, no)
        substitution_regex = "__"
        substitution_regex_to = "-"
        
        labels <- left_join(labels, tissue_names, by = "tissue") %>% left_join(., color_labels, by = "new_category")
        
      }

        labels <- labelsPrep()
    }

for_plotting <- todos %>% mutate("Name" = gsub(x = Name, pattern = "liver", replacement = "Liver")) %>% 
                                   left_join(., labels, by  = "Name") %>% arrange(new_category) 

if(args$include_factors != "ALL")
{
  source_choices = str_split(args$include_factors, pattern = ",")[[1]]
  
}else
{
  source_choices = unique(for_plotting$Source)
}
for_plotting <- for_plotting %>% filter(Source %in% source_choices) %>% mutate("z_score" = Coefficient/Coefficient_std_error)



factor_order <- paste0("F", seq(1:1000))
name_order <- unique(for_plotting$Name)
unique_categories <- unique(for_plotting$new_category)
if(FALSE)
{
  #Basic plot, nothing fancy
  ggplot(for_plotting, aes(x= factor(Source, level = factor_order), y = factor(Name, level = name_order), fill = `-log10(P)`)) + geom_tile() + scale_fill_gradient2(low="blue", high="red", mid = "white") + scale_x_discrete(guide = guide_axis(n.dodge = 1))

library(scales)
library(RColorBrewer)
#colors <- brewer.pal(length(unique_categories), "Paired")
#as of 10/18/2022, I've manually specified colors, so this category may not matter
colors <-colorRampPalette(brewer.pal(12, "Paired"))(length(unique_categories))
color_key <- data.frame("new_category" = unique_categories, "Color" = colors) #%>% rename("Tissue category for display" = `Tissue.category.for.display`)
for_plotting <- left_join(for_plotting, color_key)
color_ <- for_plotting %>% group_by(Name) %>% filter(row_number()==1) %>% select(Name, Color) %>% mutate("Name" = factor(Name, level = name_order)) %>% arrange(Name) 
}
for_plotting_strict <- for_plotting %>% 
  mutate("fdr" = p.adjust(for_plotting$Coefficient_P_value, method = "fdr"), "bonf" = p.adjust(for_plotting$Coefficient_P_value, method = "bonferroni")) %>% 
  mutate("bonfpass" = ifelse(bonf < 0.01, 1, 0), "fdrpass" = ifelse(fdr < args$fdr, 1, 0))

#If we want a tissue-specific factor or global correction:
fdr_thresh = as.character(args$fdr)
if(args$plot_type == "data_tab")
{
  for_plotting <- for_plotting  %>% mutate("global_fdr"=p.adjust(Coefficient_P_value, method="fdr")) %>%
    group_by(Source) %>% mutate("Factor_specific_fdr"=p.adjust(Coefficient_P_value, method="fdr")) %>% ungroup() %>%
    group_by(tissue, Source) %>% mutate("ts_bonfp" = p.adjust(Coefficient_P_value, method = "bonferroni")) %>% 
    slice_min(ts_bonfp, n = 1, with_ties = FALSE)  %>% ungroup() %>% group_by(Source) %>%
    mutate("factor_tissue_fdr" = p.adjust(ts_bonfp, method = "fdr")) %>% ungroup()
  write.csv(for_plotting,file=args$output,quote = FALSE, row.names = FALSE)
  quit()
}

print(args$plot_type)
if(args$plot_type == "factor_tissue_FDR" | args$plot_type == "global_tissue_FDR"){ 
  wid = 21
  ht = 14
  if(args$plot_type == "factor_tissue_FDR"){
    for_plotting_strict <- for_plotting %>% group_by(tissue, Source) %>% 
      mutate("ts_bonfp" = p.adjust(Coefficient_P_value, method = "bonferroni")) %>% 
      slice_min(ts_bonfp, n = 1, with_ties = FALSE)  %>% ungroup() %>% group_by(Source) %>%
      mutate("fdr" = p.adjust(ts_bonfp, method = "fdr")) %>% 
      ungroup() %>% mutate("fdrpass" = ifelse(fdr < args$fdr, 1, 0)) %>%
      mutate("color_code"=ifelse(fdr < 0.05, color_code,"#FFFFFF"))
    title_ = paste0("Tissue-specific enrichment,\nfactor tissue adjusted FDR < ", fdr_thresh)
    #args$plot_type = "fdr_sig"
  }
  
  if(args$plot_type == "global_tissue_FDR")
  {
    for_plotting_strict <- for_plotting %>% group_by(tissue, Source) %>% 
      mutate("ts_bonfp" = p.adjust(Coefficient_P_value, method = "bonferroni")) %>% 
      slice_min(ts_bonfp, n = 1, with_ties = FALSE) %>% ungroup() %>%
      mutate("fdr" = p.adjust(ts_bonfp, method = "fdr")) %>% 
      mutate("fdrpass" = ifelse(fdr < args$fdr, 1, 0)) %>%
      mutate("color_code"=ifelse(fdr < 0.05, color_code,"#FFFFFF"))
    title_ = paste0("Tissue-specific enrichment,\nglobal tissue adjusted FDR < ", fdr_thresh)
    
    #args$plot_type = "fdr_sig"
  }
  tissue.order <- unique((for_plotting_strict %>% arrange(new_category))$p_tissue)
  
  pass_tissues <- unique((for_plotting_strict %>% filter(fdr < args$fdr))$tissue)
  pass_cat_names <- unique((for_plotting_strict %>% arrange(color_code) %>% filter(fdr < args$fdr))$new_category)
  pass <- for_plotting_strict %>% filter(tissue %in% pass_tissues) %>% mutate("coef_out" = ifelse(fdrpass == 1, z_score, 0 ))
 save(pass, pass_tissues, pass_cat_names, file = "/scratch16/abattle4/ashton/snp_networks/presentation_figures/ashg_2023/tissue_enrich.MAY_VERSION_0.05.RData")
  p = ggplot(pass, aes(x= factor(Source, level = factor_order), y =factor(p_tissue, level = tissue.order), 
                   fill = color_code, alpha = coef_out)) + 
    geom_tile(color = "gray")  + scale_fill_identity(guide = "legend", labels =pass_cat_names) + 
    scale_y_discrete(label=function(x) abbreviate(gsub(x,pattern = "_", replacement = " "), minlength = 20)) + 
    ylab("Tissue type") + xlab("Factor Number") + ggtitle(title_) +
    guides(fill=guide_legend(title="Tissue Category")) + theme_classic(17) + labs("alpha" = "Z-score") 
  
  if(args$orientation == "horizontal")
  {
    p <- p + coord_flip() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") + 
      theme(panel.grid.major.x = element_blank())
  }
}else if(args$plot_type == "vertical")
{ 
  p = ggplot(for_plotting, aes(x= factor(Source, level = factor_order), y = factor(Name, level = name_order), fill = `-log10(P)`)) + geom_tile() + 
    scale_fill_gradient2(low="blue", high="red", mid = "white") + scale_x_discrete(guide = guide_axis(n.dodge = 1)) + theme(axis.text.y = element_text(colour=color_$Color, angle = 25)) +
    scale_y_discrete(label=function(x) abbreviate(gsub("A\\d+\\.(\\d\\d\\d.)*", "", x),minlength = 7)) + 
    ylab("Tissue type") + xlab("Factor Number")
  
} else if (args$plot_type == "horizontal") {
    wid = 15
    ht = 6 
    #print("Making horizontal plot") 
    p = ggplot(for_plotting, aes(x= factor(Source, level = factor_order), y = factor(Name, level = name_order), alpha = `-log10(P)`, fill = new_category)) +
    geom_tile() + coord_flip() + 
    theme(axis.text.x = element_text(angle = 90)) + 
    scale_y_discrete(label=function(x) abbreviate(gsub("A\\d+\\.(\\d\\d\\d.)*", "", x),minlength = 7)) + 
    ylab("Tissue type") + xlab("Factor Number") + guides(fill = guide_legend(title="Tissue Category"))

} else if (args$plot_type == "facet_grid") {
  
  p=ggplot(for_plotting, aes(x= factor(Source, level = factor_order), y = factor(Name, level = name_order), fill = `-log10(P)`)) + geom_tile() + 
    scale_fill_gradient2(low="blue", high="red", mid = "white") + ylab("Tissue type") + xlab("Factor Number") + facet_grid(~new_category) + 
    theme(axis.text.y=element_blank())
  
} else if (args$plot_type == "facet_wrap") {
  wid = 15
  ht = 7
  
  p=ggplot(for_plotting , aes(x = Name, y= factor(Source, level = factor_order), fill = Coefficient/Coefficient_std_error)) + 
    geom_tile(color = "gray")  + 
    facet_grid(~new_category, scales = "free_x", switch = "x", space = "free_x") + theme_minimal(15) + 
    scale_fill_gradient2(low="blue", high="red", mid = "white") + theme(strip.text = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Factor") + labs(fill="Z score") + theme(axis.title = element_blank(),axis.text.x=element_blank(), axis.ticks.x = element_blank())
  
  
  
}else if (args$plot_type == "bonf_sig") {
  
  #Plot just those that reach a bonferroni significance level, and then ones that reach an fdr
  p = ggplot(for_plotting_strict, aes(x= factor(Source, level = factor_order), y = factor(Name, level = name_order), fill = bonfpass)) + geom_tile(color = "gray") + 
   scale_fill_gradient2(low="blue", high="red", mid = "white") + coord_flip() + theme(axis.text.x = element_text(colour=color_$Color, angle = 90)) + 
   scale_y_discrete(label=function(x) abbreviate(gsub("A\\d+\\.(\\d\\d\\d.)*", "", x),minlength = 7)) + ylab("Tissue type") + xlab("Factor Number")
  
}else if (args$plot_type == "fdr_sig") {
    wid = 15
    ht = 8
  pass_tissues <- unique((for_plotting_strict %>% filter(fdr < args$fdr))$Name)
  fdr_thresh = as.character(args$fdr)
    if(length(pass_tissues) == 0)
    {
        print("FDR threshold too stringent. Adjusting to 0.1")
        pass_tissues <- unique((for_plotting_strict %>% filter(fdr < 0.1))$Name)
        fdr_thresh = "0.10"
        print("passing tissues")
        print(pass_tissues)
    }
    pass <- for_plotting_strict %>% filter(Name %in% pass_tissues)
    #name_order <- unique(pass$Name)
    name_order <- unique(pass$Name_source)
    if(is.null(name_order)){ name_order <- unique(pass$Name) }
    if(args$include_factors != "ALL")
    {
      #if just a subplot, do it vertically.
      p = ggplot(pass, aes(x= factor(Source, level = factor_order), y = factor(Name, level = name_order), fill = factor(new_category), alpha = fdrpass)) + 
        geom_tile(color = "gray")  + 
        scale_y_discrete(label=function(x) abbreviate(gsub(substitution_regex, substitution_regex_to, x),minlength = 25)) + 
        ylab("Tissue type") + xlab("Factor Number") + ggtitle(paste0("Tissue-specific enrichment, FDR < ", fdr_thresh)) +
        guides(fill=guide_legend(title="Tissue Category"), alpha = "none") + theme_minimal(15)
    }else
    {
      p = ggplot(pass, aes(x= factor(Source, level = factor_order), y = factor(Name, level = name_order), fill = factor(new_category), alpha = fdrpass)) + 
        geom_tile(color = "gray")  + coord_flip() + 
        scale_y_discrete(label=function(x) abbreviate(gsub(substitution_regex, substitution_regex_to, x),minlength = 25)) + 
        ylab("Tissue type") + xlab("Factor Number") + ggtitle(paste0("Tissue-specific enrichment, Global FDR < ", fdr_thresh)) +
        guides(fill=guide_legend(title="Tissue Category"), alpha = "none") + theme_classic(17) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +theme(panel.grid.major.x = element_blank())
    }

    
    #For my presentation plots
    if(FALSE)
    {
      p = ggplot(pass, aes(x= factor(Source, level = factor_order), y = factor(Name, level = name_order), fill = factor(new_category), alpha = fdrpass)) + 
        geom_tile(color = "gray") + coord_flip()  + 
        scale_y_discrete(label=function(x) abbreviate(gsub(substitution_regex, substitution_regex_to, x),minlength = 50)) + 
        ylab("Tissue type") + xlab("Factor Number") + ggtitle(paste0("Tissue-specific loading enrichment, FDR < ", fdr_thresh)) +
        guides(fill=guide_legend(title="Tissue Category"), alpha = FALSE) + theme_minimal(17) + 
        theme(axis.text.x = element_text(angle = 45, hjust =1), legend.position = "bottom")
    }

} else {
  #Leftover code- FDR at 0.1
  associations_per_factor.0.5 <- for_plotting_strict %>% group_by(Source) %>% summarize("tot" = sum(fdrpass))
  #If we loosen it up....
  for_plotting_strict <- for_plotting %>% 
    mutate("fdr" = p.adjust(for_plotting$Coefficient_P_value, method = "fdr"), "bonf" = p.adjust(for_plotting$Coefficient_P_value, method = "bonferroni")) %>% 
    mutate("bonfpass" = ifelse(bonf < 0.01, 1, 0), "fdrpass" = ifelse(fdr < 0.10, 1, 0))
  
  pass_tissues <- unique((for_plotting_strict %>% filter(fdr < 0.1))$Name)
  pass <- for_plotting_strict %>% filter(Name %in% pass_tissues)
  message("nothing left...")
  p = ggplot(pass, aes(x= factor(Source, level = factor_order), y = factor(Name, level = name_order), fill = factor(new_category), alpha = fdrpass)) + 
    geom_tile(color = "gray")  + coord_flip() + 
    theme(axis.text.x = element_text(colour=(color_ %>% filter(Name %in% pass_tissues))$Color, angle = 90)) + 
    scale_y_discrete(label=function(x) abbreviate(gsub(substitution_regex, substitution_regex_to, x),minlength = 15)) + ylab("Tissue type") + 
    xlab("Factor Number")
  
}

ggsave(args$output, plot = p, width = wid, height = ht)
save(p, file=paste0(args$output, ".RData"))
#Bonus material- compare factors and tissue enrichment?
if(FALSE)
{
  yuan <- fread("/scratch16/abattle4/ashton/snp_networks/ts_eQTLs/output/sn_spMF_K7_a10.01_l10.1/sn_spMF_FactorMatrix_K7_a10.01_l10.1.txt")
  t <- yuan %>% pivot_longer(cols = c("1","2","3","4","5","6","7"))
  ggplot(data = t, aes(x = name, y = V1, fill = value)) + geom_tile() + scale_fill_gradient2(low="blue", high="red", mid = "white")
}

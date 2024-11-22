#Plots demonstrating cohort overlap
source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/src/postMungeUtils.R")
source("/scratch16/abattle4/ashton/snp_networks/gwas_decomp_ldsc//src/plot_functions.R")


pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, cowplot, magrittr, stringdist,reshape2)
##Filll in entries in melted matrix to make matrix symmetric, if unavailable
fillInSymmetric <- function(df)
{
  cp.df <- df
  add.list <- NULL
  cp.df$all.ids <- paste0(df$phenotype_1, "_", df$phenotype_2)
  cp.df$complements <- paste0(df$phenotype_2, "_", df$phenotype_1)
  for(i in 1:length(cp.df$complements))
  {
    complement <- cp.df$complements[i]
    if(!(complement %in% cp.df$all.ids))
    {
      add.entry <- cp.df[i,]
      ref <- add.entry$desc1 
      add.entry$desc1 <- add.entry$desc2
      add.entry$desc2 <- ref
      
      ref <- add.entry$phenotype_1 
      add.entry$phenotype_1 <- add.entry$phenotype_2
      add.entry$phenotype_2 <- ref
      add.list <- rbind(add.list, add.entry)
    }
  }
  return(rbind(cp.df, add.list))
}

#Take a long data frame and turn into a correlation matrix for clustering
clusterFromLong <- function(df)
{
  updated.mat <- dcast(df, desc1 ~ desc2, value.var = "gcov_int")
  rownames.sto<- unlist(updated.mat[,1])
  updated.mat <- as.matrix(updated.mat[,-1]) %>% apply(., 1, as.numeric)
  diag(updated.mat) <- 1
  rownames(updated.mat) <- rownames.sto; colnames(updated.mat) <- rownames.sto
  order_dat <- reorder_cormat(updated.mat)
  order_dat$cormat
}

getCorrMatrix <- function(rns, cns, lookup.tab, index.list, data.df)
{
  out.mat <- matrix(NA, nrow = length(rns), ncol=length(cns))
  rownames(out.mat) <- rns; colnames(out.mat) <- cns
  for(r in rns)
  {
    for(c in cns)
    {
      if(!is.na(out.mat[r,c]))
      {
        next
      }
      print(paste("evaluating col :", c))
      row.query <- (lookup.tab %>% filter(true_name == r))$query.names
      col.query <- (lookup.tab %>% filter(true_name == c))$query.names
      if(length(row.query) == 0 | length(col.query) == 0 |length(row.query) > 1 | length(col.query)> 1)
      {
        message("Error at ", r, " ", c)
        out.mat[r,c] <- NA
      }
      else
      {
        joint <- intersect(index.list[[row.query]], index.list[[col.query]])
        choices.cols <- filter(data.df, ID %in% joint)
        row.vals <- as.matrix(choices.cols %>% select(contains(row.query))) %>% 
          apply(., 1, function(x) mean(x,na.rm=TRUE ))
        col.vals <- as.matrix(choices.cols %>% select(contains(col.query))) %>% 
          apply(., 1, function(x) mean(x,na.rm=TRUE ))
        #Alternative DATA %>% mutate(NAYSA = rowMeans(across(starts_with("STUFF"))))
        stopifnot(length(col.vals) == length(row.vals))
        out.mat[r,c] <- cor(row.vals, col.vals)
        out.mat[c,r] <- cor(row.vals, col.vals)
      }
      
    }
    print(paste("completed row :", r))
  }
  out.mat
  
}






ukbb.neale.corrs <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/Neale_lab/genetic_correlations/geno_correlations.simplified.txt")

## Data cleanup- include only the irnt version of phenotypes, not the raw ones....
drop.raws.p1 <- !grepl(ukbb.neale.corrs$p1, pattern = "_raw")
drop.raws.p2 <- !grepl(ukbb.neale.corrs$p2, pattern = "_raw")
drops <- drop.raws.p1 & drop.raws.p2
ukbb.neale.corrs.scaled <- ukbb.neale.corrs[drops,] %>% filter(!is.na(rg),abs(rg) <= 1,abs(gcov_int) <= 1, !is.na(gcov_int)) 
#rm(ukbb.neale.corrs)

#Get the data frm neal lab, including names, etc.
ukbb.names.full <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/Neale_lab/heritability_estimates/ukb31063_h2_all.02Oct2019.tsv")[,c(1:25,778)] 

## Cleanup the columns etc.
ukbb.neale.with.n <- ukbb.neale.corrs.scaled %>% mutate("phenotype"=p1, "sex" = p1_cohort) %>% 
  left_join(., ukbb.names.full %>% select(phenotype, sex,description, n, Neff), by = c("phenotype", "sex")) %>% select(-phenotype) %>% rename("Neff_p1" = Neff,"Desc_p1" = description)
#now the 2nd phenotype
ukbb.neale.with.n %<>% mutate("phenotype"=p2, "sex" = p2_cohort) %>% 
  left_join(., ukbb.names.full %>%  select(phenotype, sex, description, n, Neff), by =c("phenotype", "sex")) %>% 
  select(-phenotype) %>% 
  rename("Neff_p2" = Neff,"Desc_p2" = description)


### First plot: relationship between sample size and gcov_int estimates
bin.df <- ukbb.neale.with.n  %>% mutate("bins" = cut_interval(abs(gcov_int),n = 10)) %>% select(p1,p2,rg,gcov_int, gcov_int_se, Neff_p1, Neff_p2,Desc_p1,Desc_p2, bins)

ggplot(bin.df %>% filter(abs(gcov_int) >= 0.5), aes(x = log10(Neff_p1), y = log10(Neff_p2), color = gcov_int)) + 
  geom_point(size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") + theme_bw()+ ggtitle("High covariance studies span many sample sizes")

### For plots exploring the SE associated with these estimates and how many across the bins, please see the may2024 version of this script

### For tables looking at adjusting for multiple hypothesis testing in rg, see the may2024 version of this script, or data in /scratch16/abattle4/ashton/snp_networks/presentation_figures/lab_meeting_april2023.RData"

semi_herit_traits <- ukbb.names.full %>% filter(h2_p < 0.01)
nonherit_traits <- ukbb.names.full %>% filter(h2_p >= 0.01)
herit.only <- ukbb.neale.with.n %>% filter(p1 %in% semi_herit_traits$phenotype, p2 %in% semi_herit_traits$phenotype)
nonherit.only <- ukbb.neale.with.n %>% filter(p1 %in% nonherit_traits$phenotype, p2 %in% nonherit_traits$phenotype)
herit.plt <- ggplot(herit.only, aes(x=gcov_int)) + geom_histogram(bins = 50)  + scale_y_log10() +
  xlab(bquote("Estimate of covariance inflation"~g[covint])) + ggtitle(bquote(p[h^2]~"<0.01")) + theme_classic()+
  ylab("log10 count of UKBB study pairs")

nonherit.plt <- ggplot(nonherit.only, aes(x=(gcov_int+1e-20))) + geom_histogram(50) + scale_y_log10()  + 
  xlab(bquote("Estimate of covariance inflation"~g[covint])) + ggtitle(bquote(p[h^2]~">0.01")) + theme_classic()+
  ylab("log10 count of UKBB study pairs")
cowplot::plot_grid(plotlist=list(herit.plt, nonherit.plt))

###Comparison to published phenotypes (Tanigawa et al, 2018)

#### For the following, consider only those gcov_int which fall into bonferoni corrected 95% CI
nintyfive <- abs(qnorm(0.05/2))
ukbb.neale.corrs.scaled  %<>% mutate("upper_CI" = gcov_int + nintyfive*gcov_int_se,"lower_CI" = gcov_int - nintyfive*gcov_int_se) %>%
  mutate("CI_contains_0" = ifelse((upper_CI >= 0 & lower_CI <= 0), TRUE, FALSE)) %>% mutate("z_score_int" = gcov_int/gcov_int_se) %>% 
  mutate("pval.int" = 2*pnorm(-abs(z_score_int)))
bonf <- abs(qnorm((0.05/nrow(ukbb.neale.corrs.scaled))/2)) 

ukbb.neale.corrs.scaled  %<>% mutate("upper_CI.bonf" = gcov_int + bonf*gcov_int_se,"lower_CI.bonf" = gcov_int - bonf*gcov_int_se) %>%
  mutate("CI_contains_0.bonf" = ifelse((upper_CI.bonf >= 0 & lower_CI.bonf <= 0), TRUE, FALSE))

filtered.cov.bonf <- ukbb.neale.corrs.scaled %>% filter(!CI_contains_0.bonf)


#### Get the phenotype data for matching
ukbb.names <- ukbb.names.full %>% 
  rename("p1" = phenotype, "p1_cohort"=sex) %>% select(p1, description,p1_cohort) 
filtered.cov.bonf.desc <- filtered.cov.bonf %>% left_join(., ukbb.names, by = c("p1", "p1_cohort")) %>% select(p1,description,p1_cohort,p2,p2_cohort,gcov_int,gcov_int_se) %>% 
  rename("phenotype_1"=p1, "p1_sex" = p1_cohort, "desc1"=description) %>% rename("p1"=p2, "p1_cohort"=p2_cohort) %>% left_join(., ukbb.names,by = c("p1", "p1_cohort")) %>% 
  rename("phenotype_2"=p1, "desc2"=description, "p2_sex"=p1_cohort)


tanigwa.phenos <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/comparing_to_degas_factorgo/tanigwa_phenos.tsv",sep = ",")
#from https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11953-9/MediaObjects/41467_2019_11953_MOESM5_ESM.xlsx

joint.list <- unique(c(tolower(filtered.cov.bonf.desc$desc1),tolower(filtered.cov.bonf.desc$desc2)))
d <- expand.grid(joint.list,tolower(tanigwa.phenos$Phenotype_name)) # Distance matrix in long form
names(d) <- c("desc.name","tanigawa.name")
d$dist <- stringdist(d$desc.name,d$tanigawa.name, method="jw") # String edit distance (use your favorite function here)
better.d <- d %>%  distinct() %>% group_by(desc.name) %>% arrange(dist) %>% slice_head() %>% ungroup() %>% filter(dist < 0.103) %>% 
  filter(tanigawa.name != "vitamin a")

## Repeat this plot, but just show the ones in tanigawa pairwise:
both.in.tanigawa <- filtered.cov.bonf.desc %>% drop_na() %>% mutate('in.tanigawa' = ifelse(tolower(desc1) %in% better.d$desc.name, TRUE, FALSE)) %>%
 mutate('in.tanigawa2' = ifelse(tolower(desc2) %in% better.d$desc.name, TRUE, FALSE)) %>% filter(in.tanigawa2,in.tanigawa)
both.in.counts <- both.in.tanigawa %>% filter(abs(gcov_int) > 0.5) %>% group_by(desc1) %>% summarise("count" = n()) %>% arrange(-count)
j <- both.in.tanigawa %>% filter(in.tanigawa, in.tanigawa2)
uniqe.phenos<- unique(c(j$phenotype_1, j$phenotype_2))

### Density plots as aabove, but just looking at those in Tanigawa with 95% bonferroni estimates

ggplot(both.in.tanigawa %>% filter(phenotype_1 %in% semi_herit_traits$phenotype, phenotype_2 %in% semi_herit_traits$phenotype), aes(x=gcov_int)) + 
  geom_histogram(bins = 50)+ scale_y_log10() +
  xlab(bquote("Estimate of covariance inflation"~g[covint])) + ggtitle(bquote(p[h^2]~"<0.01")) + theme_classic() +
  ylab("log10 count of UKBB study pairs")

ggplot(both.in.tanigawa %>% filter(phenotype_1 %in% nonherit_traits$phenotype, phenotype_2 %in% nonherit_traits$phenotype), aes(x=gcov_int)) + 
  geom_histogram(bins = 50)+ scale_y_log10() +
  xlab(bquote("Estimate of covariance inflation"~g[covint])) + ggtitle(bquote(p[h^2]~">0.01")) + theme_classic() +
  ylab("log10 count of UKBB study pairs")

### Heatmap version (again, just those with a 95% bonferroni corrected estimate )
ggplot(both.in.tanigawa, aes(x = phenotype_1, y = phenotype_2, fill = gcov_int)) + geom_tile() + theme_test(16)+
  scale_fill_gradient2(low="blue",mid = "white", high="red") + 
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  ) + xlab("GWAS 1") +scale_x_discrete(position = "top")+  ylab("GWAS 2") + 
  ggtitle("Covariance inflation estimates for 292 GWAS") + labs(fill = "Estimated\noverlap\neffect")


## Merging with those in factorgo too
those.in.factor.go <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/comparing_to_degas_factorgo/factorgo.phenos.tsv")

both.in.both <- both.in.tanigawa %>% filter(desc1 %in% those.in.factor.go$description, desc2 %in% those.in.factor.go$description)

### Group by trait category
watanabe.assignments <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/comparing_to_degas_factorgo/watanabe_assignments.tsv")
both.in.watanabe <- both.in.both %>% filter(toupper(desc1) %in% toupper(watanabe.assignments$Trait),  toupper(desc2) %in% toupper(watanabe.assignments$Trait))
easy.matches <- watanabe.assignments  %>% 
  filter(toupper(Trait) %in% toupper(both.in.watanabe$desc1),  toupper(Trait) %in% toupper(both.in.watanabe$desc2)) %>% 
  select(Trait, Domain)


lookup.list  <- unique(c(both.in.both$desc1,both.in.both$desc2))
lookup.list  <- lookup.list[!(lookup.list  %in% both.in.watanabe$desc1)]
lookup.list  <- lookup.list[!(lookup.list  %in% both.in.watanabe$desc2)]
#manuall made some things..
#Make sure we have all the data from the ukbb
lookup.list  <- unique(c(both.in.both$phenotype_1,both.in.both$phenotype_2))
name.data <- ukbb.names.full %>% filter(phenotype %in% lookup.list) %>% select(phenotype, description) %>%
  mutate("p1"=phenotype, "p2"=phenotype, "desc1"=description, "desc2"=description) %>% distinct()

data.restart <- ukbb.neale.corrs %>% filter(p1 %in% lookup.list & p2 %in% lookup.list) %>% 
  left_join(.,select(name.data, p1, desc1), by="p1" ) %>% left_join(.,select(name.data, p2, desc2), by="p2") %>%
  select(p1, desc1, p2, desc2, gcov_int, gcov_int_se, rg, p) %>% 
  set_colnames(c("phenotype_1","desc1", "phenotype_2","desc2", "gcov_int","gcov_int_se","rg","p" ))

manual.assigns <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/comparing_to_degas_factorgo/phenotypes_by_category_tanigawa_factorgo.tsv")
###STOPPED here, pickup later (stopped 9/23)
##non-genetic category data for heatmaps
#c("Activities", "Social Interactions", "Environment", "Nutritional")
non.genetic <- manual.assigns %>% filter(Domain %in% c("Environment", "Nutritional"))
non.genetic.query <- filter(data.restart, desc1 %in% non.genetic$Trait, desc2 %in% non.genetic$Trait)
#Narrow it down- drop a few manually:
omits <- c("Fresh fruit intake", "Dried fruit intake", "Bread intake", "Cooked vegetable intake",
           "Average evening sound level of noise pollution", "Average 16-hour sound level of noise pollution",
           "Nitrogen dioxide air pollution; 2007", "Traffic intensity on the nearest major road", 
           "Water intake", "Cereal intake", "Nitrogen dioxide air pollution; 2006", 
           "Average 24-hour sound level of noise pollution","Nitrogen dioxide air pollution; 2005",
           "Particulate matter air pollution (pm2.5) absorbance; 2010", "Tea intake", "Coffee intake")
non.genetic.query %<>% filter(!(desc1 %in% omits), !(desc2 %in% omits))

ggplot(non.genetic.query, aes(x = phenotype_1, y = phenotype_2, fill = gcov_int)) + geom_tile() + theme_test()+
  scale_fill_gradient2(low="blue",mid = "white", high="red") + xlab("Phenotype 1") + ylab("Phenotype 2") + 
  ggtitle("Distribution of covariance inflation varies\nbetween phenotype pairs") 

non.gene
tic.query.full <- fillInSymmetric(non.genetic.query)
non.genetic.query.corr <- clusterFromLong(non.genetic.query.full)
plotCorrelationHeatmap(non.genetic.query.corr, drop.labels=TRUE, col_limit=c(-1,1.1), font_size = 16) + 
  theme(legend.position = "none")
#What is the heritability of these selected traits
herits.environmental <- filter(ukbb.names.full, phenotype %in% non.genetic.query$phenotype_1 | phenotype %in% non.genetic.query$phenotype_2 ) %>%
  filter(sex == "both_sexes")
ggplot(herits.environmental, aes(x=reorder(description, h2_liability), y= h2_liability, fill = -log10(h2_p))) + geom_bar(stat="identity") +
  coord_flip() + theme_bw(15) + xlab("Phenotype") + scale_fill_gradient(breaks=c(2,5,10,20,30), low = "grey", high="red")
  
ggplot(herits.environmental, aes(x=reorder(description, h2_p), y= -log10(h2_p), fill = h2_liability)) + geom_bar(stat="identity") +
  coord_flip() + theme_bw(15) + xlab("Phenotype") + scale_fill_gradient(low = "grey", high="red")+ geom_hline(yintercept = -log10(0.01)) +
  ggtitle("Heritability of selected phenotypes")

#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Just be more rigorous about it
low.herit <- filter(ukbb.names.full, h2_p > 0.05) %>% filter(sex == "both_sexes") %>% 
  filter(description %in% manual.assigns$Trait) %>% filter(!grepl(x = phenotype, pattern = "_raw"))
drops.low <- c("Average 24-hour sound level of noise pollution","Psoriasis", 
               "Osteoporosis", "Tinnitus severity/nuisance",
               "Age heart attack diagnosed", "Number of stillbirths", "QRS duration", "P duration",
               "Cancer year/age first occurred", "Age at death", 
               "Lifetime number of same-sex sexual partners", "Age cataract diagnosed",
               "Age diabetes diagnosed", "Mean signal-to-noise ratio (SNR), (left)",
               "Mean signal-to-noise ratio (SNR), (right)","Years of cough on most days",
               "6mm index of best keratometry results (left)",
               "Number of unenthusiastic/disinterested episodes")
low.herit %<>% filter(!(description %in% drops.low))
no.herit.query <- filter(data.restart, desc1 %in% low.herit$description, desc2 %in% low.herit$description)
no.herit.query.full <- fillInSymmetric(no.herit.query)
no.herit.query.corr <- clusterFromLong(no.herit.query.full)
plotCorrelationHeatmap(no.herit.query.corr, drop.labels=TRUE, col_limit=c(-1,1.2), font_size = 16) + 
  theme(legend.position = "none")


ggplot(low.herit, aes(x=reorder(description, h2_p), y= (h2_p), fill = h2_liability)) + geom_bar(stat="identity") +
  coord_flip() + theme_bw(15) + xlab("Phenotype") + scale_fill_gradient(low = "grey", high="red")+ geom_hline(yintercept = 0.05) +
  ggtitle("Heritability of selected phenotypes")


#Look just at metabolic traits?
metabolic<- manual.assigns %>% filter(Domain %in% c("Metabolic"))
metabolic.query <- filter(data.restart, desc1 %in% metabolic$Trait, desc2 %in% metabolic$Trait)
#Raw
ggplot(metabolic.query, aes(x = phenotype_1, y = phenotype_2, fill = gcov_int)) + geom_tile() + theme_test()+
  scale_fill_gradient2(low="blue",mid = "white", high="red") + xlab("Phenotype 1") + ylab("Phenotype 2") + ggtitle("Distribution of covariance inflation varies\nbetween phenotype pairs") 

metabolic.query.full <- fillInSymmetric(metabolic.query)
metabolic.query.corr <- clusterFromLong(metabolic.query.full)
plotCorrelationHeatmap(metabolic.query.corr, drop.labels=TRUE, col_limit=c(-1,1.2)) + labs(fill="test")

#this is too much- just pick a few...
#Now, for some of the other traits:
manual.choice <- manual.assigns %>% filter(Domain %in% c("Metabolic","Cardiovascular", "Ophthalmological"))
choices <- c("3mm strong meridian (left)","6mm strong meridian (left)",
             "6mm strong meridian (right)","3mm strong meridian (right)","Body mass index (BMI)", "Basal metabolic rate",
             "Whole body fat mass", "Leg fat mass (left)", "Weight", "Whole body water mass", "Systolic blood pressure, manual reading", 
             "Pulse rate, automated reading", "Diastolic blood pressure, manual reading","Ventricular rate")
manual.query <- filter(data.restart, desc1 %in% choices, desc2 %in% choices)
#There are some redundant phenotypes, either rename or drop
manual.query %<>% filter(!(phenotype_1 %in% c("23104_irnt","23098_irnt")),
                         !(phenotype_2 %in% c("23104_irnt","23098_irnt")))

ggplot(manual.query, aes(x = desc1, y = desc2, fill = gcov_int)) + geom_tile() + theme_test()+
  scale_fill_gradient2(low="blue",mid = "white", high="red") + xlab("Phenotype 1") + ylab("Phenotype 2") + ggtitle("Distribution of covariance inflation varies\nbetween phenotype pairs") 

manual.query.full <- fillInSymmetric(manual.query)
manual.query.corr <- clusterFromLong(manual.query.full)
plotCorrelationHeatmap(manual.query.corr, drop.labels=TRUE, col_limit=c(-1,1.2), font_size = 16) + 
  theme(legend.position = "right")+  labs(fill=bquote(g[cov_int])) 


#updated.mat <- dcast(tops, desc1 ~ desc2, value.var = "gcov_int")
#rownames(updated.mat) <- updated.mat[,1]
#updated.mat <- as.matrix(updated.mat[,-1])
#order_dat <- reorder_cormat(updated.mat)
#cormat <- order_dat$cormat






#########
#LOOKING AT overlap to add to the plot
ukbb.phenos.part1 <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/comparing_to_degas_factorgo/gleaner_plot_part_1_participant.tsv")[,-105]
ukbb.phenos.part2 <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/comparing_to_degas_factorgo/gleaner_2_participant.tsv")
ukbb.phenos.part3 <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/comparing_to_degas_factorgo/run_3_keratometry_participant.tsv")
colnames.correct <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/comparing_to_degas_factorgo/gleaner_plot_1_ids.tsv",header = FALSE, sep = "\t")
colnames(ukbb.phenos.part1) <- colnames.correct$V1[-105]
cols.tab.2 <- c("ID","Whole body water mass | Instance 0","Whole body water mass | Instance 1",
                "Whole body water mass | Instance 2","Whole body water mass | Instance 3","Basal metabolic rate | Instance 0",
                "Basal metabolic rate | Instance 1","Basal metabolic rate | Instance 2","Basal metabolic rate | Instance 3",
                "Systolic blood pressure manual reading | Instance 0 | Array 0","Systolic blood pressure manual reading | Instance 0 | Array 1",
                "Systolic blood pressure manual reading | Instance 1 | Array 0","Systolic blood pressure manual reading | Instance 1 | Array 1",
                "Systolic blood pressure manual reading | Instance 2 | Array 0","Systolic blood pressure manual reading | Instance 3 | Array 0",
                "Systolic blood pressure manual reading | Instance 2 | Array 1","Systolic blood pressure manual reading | Instance 3 | Array 1",
                "Diastolic blood pressure manual reading | Instance 0 | Array 0","Diastolic blood pressure manual reading | Instance 0 | Array 1",
                "Diastolic blood pressure manual reading | Instance 1 | Array 0","Diastolic blood pressure manual reading | Instance 1 | Array 1",
                "Diastolic blood pressure manual reading | Instance 2 | Array 0","Diastolic blood pressure manual reading | Instance 2 | Array 1",
                "Diastolic blood pressure manual reading | Instance 3 | Array 0","Diastolic blood pressure manual reading | Instance 3 | Array 1")
colnames(ukbb.phenos.part2) <- cols.tab.2
stopifnot(all(ukbb.phenos.part1$ID == ukbb.phenos.part2$ID))
stopifnot(all(ukbb.phenos.part3$`Participant ID` == ukbb.phenos.part2$ID))
#Remove columns with word "angle" and a few other errors
ukbb.phenos.joined <- cbind(ukbb.phenos.part1, ukbb.phenos.part2[,-1], ukbb.phenos.part3[,-1]) %>% select(-contains('angle')) %>% 
  select(-contains("6mm index of best keratometry results (left)")) %>%
  select(-contains("3mm index of best keratometry results (left)"))

  

#SANDBOX
indices <- c(1,which(grepl(pattern = "Ventricular", x= colnames.correct$V1)))
v.rate <- (ukbb.phenos.part1[,..indices] %>% mutate("valid_entry" = ifelse(is.na(`Ventricular rate | Instance 2`) & is.na(`Ventricular rate | Instance 3`), 0,1)) %>%
  filter(valid_entry == 1))$ID

indices <- which(grepl(pattern = "Pulse", x= colnames.correct$V1))
pulse.automated <- ukbb.phenos.part1$ID[which(as.matrix(ukbb.phenos.part1[,..indices]) %>% apply(., 1, function(x) sum(!is.na(x))) > 0)]
sum(v.rate %in% pulse.automated)/length(pulse.automated)
sum(pulse.automated %in% v.rate)/length(v.rate)
sum(pulse.automated %in% v.rate)/ (sqrt(length(pulse.automated)) * sqrt(length(v.rate)))
#Beautious and gloriful. This is what we'd expect
#so for each group, get the valid ids.
queries <- unique(gsub(x=colnames(ukbb.phenos.joined), pattern = " \\| Instance \\d", replacement = "") %>% 
  gsub(x=., pattern = " \\| Array \\d", replacement = ""))[-1]
queries

index.list <- list()
for(q in queries)
{
  
  indices <- which(grepl(pattern = q, x= colnames(ukbb.phenos.joined),fixed=TRUE))
  print(colnames(ukbb.phenos.joined)[indices])
  index.list[[q]] <- ukbb.phenos.joined$ID[which(as.matrix(ukbb.phenos.joined[,..indices]) %>% apply(., 1, function(x) sum(!is.na(x))) > 0)]
}


correct <- c(colnames(no.herit.query.corr),colnames(manual.query.corr))
lookup.tab <- data.frame("query.names" = queries, "true_name" = correct[stringdist::amatch(queries, correct, maxDist = Inf)]) %>% 
  mutate("exact_same" = ifelse(query.names==true_name, TRUE, FALSE))
stopifnot(length(unique(lookup.tab$true_name)) == length(lookup.tab$true_name))

out.mat <- no.herit.query.corr
getCountsMatrix <- function(rns, cns, lookup.tab, index.list)
{
  out.mat <- matrix(NA, nrow = length(rns), ncol=length(cns))
  rownames(out.mat) <- rns; colnames(out.mat) <- cns
  for(r in rns)
  {
    for(c in cns)
    {
      row.query <- (lookup.tab %>% filter(true_name == r))$query.names
      col.query <- (lookup.tab %>% filter(true_name == c))$query.names
      if(length(row.query) == 0 | length(col.query) == 0 |length(row.query) > 1 | length(col.query)> 1)
      {
        message("Error at ", r, " ", c)
        out.mat[r,c] <- NA
      }
      else
      {
        overall <- sum(index.list[[row.query]] %in% index.list[[col.query]])
        stopifnot(overall == sum(index.list[[col.query]] %in% index.list[[row.query]]))
        out.mat[r,c] <- overall/(sqrt(length(index.list[[row.query]])) * sqrt(length(index.list[[col.query]])))
      }
      
    }
  }
  out.mat
  
}


overlap.no.herit <- round(getCountsMatrix(rownames(no.herit.query.corr),colnames(no.herit.query.corr),lookup.tab,index.list), digits = 2)
corr.no.herit <- round(getCorrMatrix(rownames(no.herit.query.corr),colnames(no.herit.query.corr),lookup.tab,index.list, ukbb.phenos.joined), digits = 2)

#combine into a matrix 
overlap.and.corr <- overlap.no.herit
overlap.and.corr[lower.tri(overlap.and.corr)] <- corr.no.herit[lower.tri(corr.no.herit)]
diag(overlap.and.corr) <- NA
plotCorrelationHeatmap(no.herit.query.corr, drop.labels=TRUE, col_limit=c(-1,1.2), font_size = 16,show.nums = TRUE,print.nums = overlap.and.corr) + 
  theme(legend.position = "right")+  labs(fill=bquote(g[cov_int])) 

#this looks good! 3 quick spot checks pass

#Now the main traits
overlap.selected.traits <- round(getCountsMatrix(rownames(manual.query.corr),colnames(manual.query.corr),lookup.tab,index.list), digits =2)
corr.selected <- round(getCorrMatrix(rownames(manual.query.corr),colnames(manual.query.corr),lookup.tab,index.list, ukbb.phenos.joined), digits = 2)
#combine into a single matrix
o.and.c.selected <- overlap.selected.traits
o.and.c.selected[lower.tri(o.and.c.selected)] <- corr.selected[lower.tri(corr.selected)]
diag(o.and.c.selected) <- NA
plotCorrelationHeatmap(manual.query.corr, drop.labels=TRUE, col_limit=c(-1,1.2), font_size = 16,show.nums = TRUE,print.nums = o.and.c.selected) + 
  theme(legend.position = "right")+  labs(fill=bquote(g[cov_int])) 

save(manual.query.corr, no.herit.query.corr,o.and.c.selected,overlap.and.corr,overlap.selected.traits,overlap.no.herit,corr.no.herit,corr.selected
     file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/ukbb_overlap_dat.RData")

#load("/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/ukbb_overlap_dat.RData")
#Now make those plots..
overlap.no.herit.half <- overlap.no.herit
overlap.no.herit.half[upper.tri(overlap.no.herit.half)] <- NA
noherit.corr <- plotCorrelationHeatmap(no.herit.query.corr, drop.labels=TRUE, col_limit=c(-1,1.2), font_size = 8,show.nums = TRUE,print.nums = overlap.no.herit.half) + 
  theme(legend.position = "right")+  labs(fill=bquote(g[cov_int]))  + theme(axis.text = element_blank()) + 
  annotate(geom="text",label = bquote(h^2 ~"= 0"),y = 12,x=4, size= 17/.pt)


overlap.selected.traits.half <- overlap.selected.traits
overlap.selected.traits.half[upper.tri(overlap.selected.traits.half)] <- NA
herit.corr <- plotCorrelationHeatmap(manual.query.corr, drop.labels=TRUE, col_limit=c(-1,1.2), font_size = 8,show.nums = TRUE,print.nums = overlap.selected.traits.half) + 
  theme(legend.position = "right")+  labs(fill=bquote(g[cov_int]))   + theme(axis.text = element_blank()) + 
  annotate(geom="text",label = bquote(h^2 ~"> 0"),y = 12,x=4, size= 17/.pt)



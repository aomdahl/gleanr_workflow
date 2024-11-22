#First the heatmap:
pacman::p_load(data.table, dplyr, magrittr, tidyr)

#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/R/read_in_tools.R")
#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/R/data_decorrelation.R")
library(gleanr)
library(GGally)
get_upper_tri <- function(cormat){ #From stakckoverlfow
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
} 

#Get the finngen stuff
fg.args<- list();
fg.args$covar_matrix <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark_2/summary_data/gcov_int.tab.csv"
fg.args$sample_sd <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/sample_sd_report.tsv"
fg.args$block_covar <- 0.2
fg.args$WLgamma <- 0
namin <- colnames(fread(fg.args$covar_matrix))
covar.dat.fg <- SampleOverlapCovarHandler(fg.args, namin, matrix(0,nrow=1655,ncol=42))
#Get the UKBB stuff
uk.args<- list();
uk.args$covar_matrix <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/ukbb_benchmark_2/summary_data/gcov_int.tab.csv"
uk.args$sample_sd <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/ukbb_benchmark_2/sample_sd_report.tsv"
uk.args$block_covar <- 0.2
uk.args$WLgamma <- 0
namin <- colnames(fread(fg.args$covar_matrix))
covar.dat.uk <- SampleOverlapCovarHandler(uk.args, namin, matrix(0,nrow=1655,ncol=42))




#If we want to just visualize the one that goes as input
diff.mat <- covar.dat.fg$C-covar.dat.uk$C
keep.rows <- which(apply(diff.mat,1, function(x) sum(x != 0) > 0))
nice.names <- fread("/scratch16/abattle4/ashton/snp_networks/scratch/finngen_v_ukbb/trait_mapping_pretty.csv",header = FALSE)
colnames(diff.mat) <- gsub(colnames(diff.mat),pattern = "_", replacement=" ")
rownames(diff.mat) <- gsub(colnames(diff.mat),pattern = "_", replacement=" ")

aligned.names <- left_join(data.frame("V1"=colnames(diff.mat)), nice.names,by="V1")
colnames(diff.mat) <- aligned.names$V2; rownames(diff.mat) <- aligned.names$V2
#This thing isn't what I wnat anymore.
#ggcorr(diff.mat[keep.rows,keep.rows],cor_matrix = diff.mat[keep.rows,keep.rows], low = "blue", mid = "white", high = "red", hjust = 1, size = 5,layout.exp = 10)+
#  theme_classic()

# Convert to a data frame for plotting
cor_df <- reshape2::melt(diff.mat[keep.rows,keep.rows])

# Rename the columns for clarity
colnames(cor_df) <- c("Var1", "Var2", "Correlation")

# Keep only lower triangular values (including diagonal)
cor_df <- cor_df %>%
  filter(as.numeric(Var2) >= as.numeric(Var1))

#better trait labels:

# Create the heatmap plot

lines <- c( "first" ,
            "second")

#bquote("LDSR" ~ g[covint]*" \n(after shrinkage)")
library(scales)
tri.heatmap <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "grey") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name =bquote(atop("LDSR" ~ g[covint], "(after shrinkage)"))) +
  theme_classic(15) +
  theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.line.x = element_blank(), axis.title.y=element_blank(),axis.ticks.y = element_blank()) +
  coord_fixed() +  # Ensure square tiles
  theme(legend.position = c(0.65, 0.2)) + theme(legend.title = element_text(hjust = 1,size=11,vjust=0.5),  # Right align the legend title text
                                               legend.text = element_text(hjust = 1),
                                               legend.key.height = unit(0.65, "cm"),  # Adjust height of legend keys
                                               legend.key.width = unit(0.5, "cm"),   # Adjust width of legend keys
                                               legend.box.spacing = unit(1, "lines"),
                                               legend.key= element_rect(colour = "transparent"),
                                             legend.background = element_blank()) +  # Adjust space between legend and plot
  guides(fill = guide_colorbar(title.position = "left")) +
 theme(axis.text.y =  element_text(size = 16), legend.title = element_text(size=16))
#Just the heatmap, taller
ggsave(plot = tri.heatmap, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig2_panelC_only.png",
       width = 8,height=10,units="in")
ggsave(plot=tri.heatmap, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig2_panelC_only.svg",
       width = 8,height=10,units="in")


##Supplemental figure S3- raw differences
tot.fg <-as.matrix(fread(fg.args$covar_matrix))
tot.uk <-as.matrix(fread(uk.args$covar_matrix))
stopifnot(all(colnames(tot.fg) == colnames(tot.uk)))
tot.diff <- tot.fg - tot.uk
#Update names to be clean
colnames(tot.diff) <- gsub(colnames(tot.diff),pattern = "_", replacement=" ")
rownames(tot.diff) <- gsub(colnames(tot.diff),pattern = "_", replacement=" ")
ggcorr(tot.diff, low = "blue", mid = "white", high = "red", hjust = 1, size = 5, layout.exp = 10, cor_matrix=tot.diff)


# Convert to a data frame for plotting
diff_df <- reshape2::melt(tot.diff)
colnames(diff_df) <- c("Var1", "Var2", "Correlation_diff")
# Keep only lower triangular values (including diagonal)
diff_df <- diff_df %>%
  filter(as.numeric(Var2) >= as.numeric(Var1))
lines <- c( "first" ,
            "second")
heatmap.diff.raw <- ggplot(diff_df, aes(x = Var1, y = Var2, fill = Correlation_diff)) +
  geom_tile(color = "grey") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name =bquote(atop("Difference in", g[covint]))) +
  theme_classic(15) +
  theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.line.x = element_blank(), axis.title.y=element_blank(),axis.ticks.y = element_blank()) +
  coord_fixed() +  # Ensure square tiles
  theme(legend.position = c(0.65, 0.2)) + theme(legend.title = element_text(hjust = 1,size=5,vjust=0.5),  # Right align the legend title text
                                                legend.text = element_text(hjust = 1),
                                                legend.key.height = unit(0.65, "cm"),  # Adjust height of legend keys
                                                legend.key.width = unit(0.5, "cm"),   # Adjust width of legend keys
                                                legend.box.spacing = unit(1, "lines"),
                                                legend.key= element_rect(colour = "transparent"),
                                                legend.background = element_blank()) +  # Adjust space between legend and plot
  guides(fill = guide_colorbar(title.position = "left")) +
  theme(axis.text.y =  element_text(size = 16), legend.title = element_text(size=16))

ggsave(plot = heatmap.diff.raw, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/suppl3_heatmap_diff.png",
       width = 13,height=15,units="in")
ggsave(plot=heatmap.diff.raw, filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/suppl3_heatmap_diff.svg",
       width = 10,height=13,units="in")

######################################################

######################################################
### Looking directly at correlation:
source("/scratch16/abattle4/ashton/snp_networks/scratch/finngen_v_ukbb/src/matrix_comparison_utils.R")
fg.std_block.path <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/covar_influence/"
uk.std_block.path <- "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/conservative_1e-5/covar_influence/"
std_block.results <- loadGLEANERResultsLists(fg.std_block.path, uk.std_block.path)

#Compile and scale all the data for easy access
scaled.case <- list("finngen"=list(), 'ukbb' = list())
scal  <- std_block.results$finngen$`BIC-sklearn_eBIC_K-41`
scaled.case$finngen$BIC_adj_K41$V <- apply(scal$V, 2, scale, center=FALSE); scaled.case$finngen$BIC_adj_K41$U <- apply(scal$U, 2, scale, center=FALSE); 

scaled.case$ukbb$BIC_adj_K41 <- std_block.results$ukbb$`BIC-sklearn_eBIC_K-41`
scaled.case$ukbb$BIC_adj_K41$V <- apply(scaled.case$ukbb$BIC_adj_K41$V, 2, scale, center=FALSE); scaled.case$ukbb$BIC_adj_K41$U <- apply(scaled.case$ukbb$BIC_adj_K41$U, 2, scale, center=FALSE); 


scaled.case$finngen$BIC_no_adj_K41 <- std_block.results$finngen$`BIC-sklearn_eBIC_K-41_no_covar`
scaled.case$finngen$BIC_no_adj_K41$V <- apply(scaled.case$finngen$BIC_no_adj_K41$V,2, scale, center=FALSE);  scaled.case$finngen$BIC_no_adj_K41$U <- apply(scaled.case$finngen$BIC_no_adj_K41$U, 2, scale, center=FALSE); 


scaled.case$ukbb$BIC_no_adj_K41 <- std_block.results$ukbb$`BIC-sklearn_eBIC_K-41_no_covar`
scaled.case$ukbb$BIC_no_adj_K41$V <- apply(scaled.case$ukbb$BIC_no_adj_K41$V,2, scale, center=FALSE);  scaled.case$ukbb$BIC_no_adj_K41$U <- apply(scaled.case$ukbb$BIC_no_adj_K41$U, 2, scale, center=FALSE); 

#Numbers in text
#Dropped singletons, finngen
singletons.fg.noadj <- apply(scaled.case$finngen$BIC_no_adj_K41$V, 2, function(x) which(x!=0))
singleton.factors.fg.noadj <- scaled.case$finngen$BIC_no_adj_K41$trait.names[unlist(singletons.fg.noadj[which(sapply(singletons.fg.noadj,length) == 1)])]

singletons.fg.adj <- apply(scaled.case$finngen$BIC_adj_K41$V, 2, function(x) which(x!=0))
singleton.factors.fg.adj <- scaled.case$finngen$BIC_adj_K41$trait.names[unlist(singletons.fg.adj[which(sapply(singletons.fg.adj,length) == 1)])]
sum(!(singleton.factors.fg.noadj %in% singleton.factors.fg.adj))

#Dropped singletons, UKBB
singletons.uk.noadj <- apply(scaled.case$ukbb$BIC_no_adj_K41$V, 2, function(x) which(x!=0))
singleton.factors.uk.noadj <- scaled.case$ukbb$BIC_no_adj_K41$trait.names[unlist(singletons.uk.noadj[which(sapply(singletons.uk.noadj,length) == 1)])]

singletons.uk.adj <- apply(scaled.case$ukbb$BIC_adj_K41$V, 2, function(x) which(x!=0))
singleton.factors.uk.adj <- scaled.case$ukbb$BIC_adj_K41$trait.names[unlist(singletons.uk.adj[which(sapply(singletons.uk.adj,length) == 1)])]
sum(!(singleton.factors.uk.noadj %in% singleton.factors.uk.adj))

#Get the correlation across the versions

#Adjusted, standard gleaner

v.dat <- prepMatricesForAnalysis(scaled.case$finngen$BIC_adj_K41$V,scaled.case$ukbb$BIC_adj_K41$V); 
u.dat <- prepMatricesForAnalysis(scaled.case$finngen$BIC_adj_K41$U,scaled.case$ukbb$BIC_adj_K41$U)
stopifnot(!v.dat$swap)
fg.v <- v.dat$lead; uk.v <- v.dat$second
fg.u <- u.dat$lead; uk.u <- u.dat$second
adj.version <- greedyMaxPairedCor(fg.u,fg.v, uk.u,uk.v)

uk.u.new = matrixSignsProduct(uk.u[,adj.version$order.pred], adj.version$signs)

uk.v.new = matrixSignsProduct(uk.v[,adj.version$order.pred], adj.version$signs)
adj.v <- cor.test(y=c(uk.v.new), x=c(fg.v))
adj.u <- cor.test(y=c(uk.u.new), x=c(fg.u))

#Undjusted, non-standard gleaner
v.dat <- prepMatricesForAnalysis(scaled.case$finngen$BIC_no_adj_K41$V,scaled.case$ukbb$BIC_no_adj_K41$V); 
u.dat <- prepMatricesForAnalysis(scaled.case$finngen$BIC_no_adj_K41$U,scaled.case$ukbb$BIC_no_adj_K41$U)
stopifnot(!v.dat$swap)
fg.v <- v.dat$lead; uk.v <- v.dat$second
fg.u <- u.dat$lead; uk.u <- u.dat$second
no_adj.version <- greedyMaxPairedCor(fg.u,fg.v, uk.u,uk.v)

uk.u.non = matrixSignsProduct(uk.u[,no_adj.version$order.pred], no_adj.version$signs)

uk.v.non = matrixSignsProduct(uk.v[,no_adj.version$order.pred], no_adj.version$signs)
nonadj.v <- cor.test(y=c(uk.v.non), x=c(fg.v))
nonadj.u <- cor.test(y=c(uk.u.non), x=c(fg.u))

#### are these statistically significant differences?
#Note that the best test would be a little more nuanced, since these groups really aren't entirely independent. But this is okay for now.
library(cocor)
cocor.indep.groups(r1.jk=adj.u$estimate, r2.hm=nonadj.u$estimate, n1=length(c(uk.u.new)), n2=length(c(uk.u.non)), alternative="greater", alpha=0.05, conf.level=0.95, null.value=0)
cocor.indep.groups(r1.jk=adj.v$estimate, r2.hm=nonadj.v$estimate, n1=length(c(uk.v.new)), n2=length(c(uk.v.non)), alternative="greater", alpha=0.05, conf.level=0.95, null.value=0)

#### Make the plot
calc.cor <- data.frame("type" = c("GLEANER", "GLEANER","GLEANER\n(unadjusted)","GLEANER\n(unadjusted)"), "Matrix"=c("U","V", "U","V"),
                       "estimate"=c(adj.u$estimate,adj.v$estimate,nonadj.u$estimate, nonadj.v$estimate),
                       "upper_ci"=c(adj.u$conf.int[2],adj.v$conf.int[2],nonadj.u$conf.int[2], nonadj.v$conf.int[2]),
                       "lower_ci"=c(adj.u$conf.int[1],adj.v$conf.int[1],nonadj.u$conf.int[1], nonadj.v$conf.int[1])
)

library(RColorBrewer)
myColors <- brewer.pal(3,"Oranges")[c(3,2)] 
names(myColors) <- levels(calc.cor$type)
colScale <- scale_colour_manual(name = "Method",values = myColors)
fillScale <- scale_fill_manual(name = "Method",values = myColors)

#add the lines
lines.add <- calc.cor %>% filter(type=="GLEANER") %>% select(Matrix, estimate)

orange.bars <- ggplot(calc.cor, aes(x=type,y=estimate, fill=type)) +geom_bar(stat='identity',width = 0.45) + facet_wrap(~Matrix) +
  geom_errorbar( aes(x=type, ymin=lower_ci, ymax=upper_ci), width=0.3, colour="black", alpha=0.9) + 
  theme_bw(15) + theme(strip.background =element_rect(fill="white"), axis.title.y=element_blank()) +
  geom_hline(data=lines.add, aes(yintercept=estimate), lty="dashed",col="grey")+
  ylab("Pearson correlation")  + fillScale + xlab("Covariance adjustment method") + 
  coord_flip(ylim=c(0.5,0.9))
  

save(tri.heatmap,orange.bars, file="/scratch16/abattle4/ashton/snp_networks/presentation_figures/overlap_former_fig3.Rdata")
bottom.plot <- cowplot::plot_grid(plotlist=list(tri.heatmap, orange.bars + theme(legend.position = "none")),rel_widths = c(1,0.95))
ggsave(filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig2_panelC-D.svg",
       width = 15,height=4,units="in")
ggsave(filename="/scratch16/abattle4/ashton/snp_networks/presentation_figures/gleaner_manuscript_draft/figures/fig2_panelC-D.png",
       width = 15,height=4,units="in")


###### Stats in text:
###Average difference per trait-trait pair
fg.full.mat <- as.matrix(fread(fg.args$covar_matrix))
uk.full.mat <- as.matrix(fread(uk.args$covar_matrix))
stopifnot(all(colnames(fg.full.mat) == colnames(uk.full.mat)))
mean(abs(c(fg.full.mat- uk.full.mat)))
########And some numbers in there to include:
unadj.fg.x <- std_block.results$finngen$`BIC-sklearn_eBIC_K-41_no_covar`$U %*% t(std_block.results$finngen$`BIC-sklearn_eBIC_K-41_no_covar`$V)
adj.fg.x <- std_block.results$finngen$`BIC-sklearn_eBIC_K-41`$U %*% t(std_block.results$finngen$`BIC-sklearn_eBIC_K-41`$V)

unadj.uk.x <- std_block.results$ukbb$`BIC-sklearn_eBIC_K-41_no_covar`$U %*% t(std_block.results$ukbb$`BIC-sklearn_eBIC_K-41_no_covar`$V)
adj.uk.x <- std_block.results$ukbb$`BIC-sklearn_eBIC_K-41`$U %*% t(std_block.results$ukbb$`BIC-sklearn_eBIC_K-41`$V)


unad.rrmse <- rrmse(unadj.uk.x, unadj.fg.x)
adj.rrmse <- rrmse(adj.uk.x, adj.fg.x)

(adj.rrmse-unad.rrmse)/unad.rrmse
#  4% improvement
#Do I need to jump in on those points?

#change in singletons:
singleton.changes.fg <-
  sum(apply(special.case$finngen$BIC_adj_K41$V, 2, function(x) sum(x!=0) == 1)) - 
  sum(apply(special.case$finngen$BIC_no_adj_K41$V, 2, function(x) sum(x!=0) == 1))


singleton.changes.uk <-
  sum(apply(special.case$ukbb$BIC_adj_K41$V, 2, function(x) sum(x!=0) == 1)) - 
  sum(apply(special.case$ukbb$BIC_no_adj_K41$V, 2, function(x) sum(x!=0) == 1))

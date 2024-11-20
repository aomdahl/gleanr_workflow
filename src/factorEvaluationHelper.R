

#Function to plot teh relationship between a single factor value in U and the z scores corresponding
#meant for a singleton, but could be used otherwise
#u.ref: U mnatrix with snp labels
#z.scores.all: zscores for each trait, named, with snp labels
#trait to evaluate
#factor it goes with
singletonLookup <- function(u.ref, z.scores.all, trait, factor)
{
  trait.lookup <- u.ref %>% arrange(-abs(!!sym(factor))) %>% select(SNPs, !!sym(factor)) %>% mutate("rank"=row_number())
  rank.joined <- left_join(trait.lookup,ranked.z.scores, by="SNPs")
  cortest <- cor.test(rank.joined$rank, rank.joined$rank_z, method="spearman")
  ggplot(rank.joined, aes(x=.data[[trait]], y=.data[[factor]])) + geom_point() + theme_classic() + xlab(paste0(trait, " z-score")) + 
    annotate(geom="text", label=paste0("Spearman:", round(cortest$estimate,digits=2)),x=-7,y=min(rank.joined[[factor]]))
}

#Plot barplots for all factors loaded on a given trait.
traitBarplots <- function(trait, scaled.V,xlab=FALSE,size=15,factor_list=c())
{
  barplots = list()
  whr.factors <- names(which(unlist((scaled.V %>% dplyr::filter(Traits==trait))[,-1]) !=0))
  if(length(factor_list) > 0)
  {
    whr.factors = paste0("V", factor_list)
  }
  #barplot its inclusion in each, scaled
  
  sub.v <- scaled.V %>% select(Traits, whr.factors) %>% mutate("of_choice"=ifelse(Traits==trait, "red","black"))
  for(f in whr.factors)
  {
    barplots[[f]] <- ggplot(sub.v, aes(x=Traits, y=.data[[f]], fill=of_choice)) + 
                        geom_bar(stat="identity")+ scale_fill_manual(values = c("grey60", "coral")) + theme_classic(size) +
            theme(legend.position = "none", axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) + geom_hline(yintercept = 0,color="black") + ylab(gsub("V",replacement = "F",f))
  }
  if(!xlab)
  {
    for(f_i in 1:(length(whr.factors)-1))
    {
      f <- whr.factors[f_i]
      barplots[[f]] <- dropXLab(barplots[[f]])
    }
  }
  barplots
}


plotByFactorNumber <- function(factor.vals, factor.name,  flip_coord = TRUE, all_nonzero = FALSE)
{
  if(all_nonzero)
  {
    p <- ggplot(factor.vals %>% dplyr::filter(!!sym(factor.name) !=0), aes(x = reorder(Traits,!!sym(factor.name)), y= !!sym(factor.name))) + geom_bar(stat= "identity") + xlab("Traits") + theme_bw()
    
  }else
  {
    p <- ggplot(factor.vals %>% dplyr::filter(abs(!!sym(factor.name)) > mean(abs(!!sym(factor.name)))), aes(x = reorder(Traits,!!sym(factor.name)), y= !!sym(factor.name))) + geom_bar(stat= "identity") + xlab("Traits") + theme_bw()
    
  }
  if(flip_coord)
  {
    return(p + coord_flip())
  }
  return(p)
}
plotSNPsByFactorNumber <- function(factor.vals, factor.name, show.n=10)
{
  ggplot(factor.vals %>% slice_max(abs(!!sym(factor.name)), n=show.n), aes(x = reorder(SNPs,!!sym(factor.name)), y= !!sym(factor.name))) + geom_bar(stat= "identity") + xlab("SNP") + coord_flip() + theme_bw()
}

######Look at the changes in sparisty over time
######
#first phase
plotSparsityChanges <- function(bic.dat, ret)
{
  sparsity.fitting <- rbind(
    data.frame("iter" = 1:length(bic.dat$rec.dat$U_sparsities), "U_sparsity"=bic.dat$rec.dat$U_sparsities,
               "V_sparsity"=bic.dat$rec.dat$V_sparsities, "step"="model_selection step"),
    data.frame("iter" = 1:length(ret$U_sparsities), "U_sparsity"=ret$U_sparsities,
               "V_sparsity"=ret$V_sparsities, "step"="fitting step")) %>%
    pivot_longer(cols =c("U_sparsity", "V_sparsity"),names_to = "matrix" )
  sparsity.fitting$step = factor(sparsity.fitting$step, levels = c("model_selection step","fitting step"))
  
  ggplot2::ggplot(sparsity.fitting, aes(x=iter, y=value, color=matrix)) + geom_point() + 
    facet_wrap(~step) + theme_bw() + xlab("Iteration") + ylab("Proportion of 0s")
#scales = "free"
  }

plotNormChanges <- function(bic.dat, ret)
{
  
  fitting.step.sizes <- data.frame("iteration" = 1:length(ret$Vs), 
                                   "v_norm"=sapply(ret$Vs, function(x) norm(x,"F")),
                                   "u_norm"=sapply(ret$Us, function(x) norm(x,"F")), "step"='fitting step')
  
  
  pre.fitting.sizes <- data.frame("iteration" = 1:length(bic.dat$convergence.options), 
                                  "v_norm"=sapply(bic.dat$convergence.options, function(x) x$V.norm),"u_norm"=sapply(bic.dat$convergence.options, function(x) x$U.norm), "step"='model_selection step')
  
  
  joined <- rbind(fitting.step.sizes, pre.fitting.sizes) %>% pivot_longer(cols=c("v_norm", "u_norm")) %>%
    mutate("step"=factor(step, levels = c("model_selection step","fitting step")))
  ggplot(joined, aes(x=iteration,y=log10(value),color=name)) + geom_point() + facet_wrap(~step) + theme_bw() + ylab("log(matrix norm)")
}

##Lambda max changes from run to run, right

## Estimate effective number of SNPs
estMe <- function(U)
{
  n=nrow(U)
  apply(U,2, function(x) (3*n)/moments::kurtosis(x))
}
  
calcMt <- function(U)
{
  apply(U,2, function(x) sum(x!=0))
}

#fin=factors of interest
getPolygenicity <-function(U)
{
  data.frame("Factors"=paste0("U",1:ncol(U)), "Me"=estMe(U),"Mt"=calcMt(U))
}

#plotPolygenicity(ret$U, paste0("U",shbg.factors),type = "prop")
plotPolygenicity <- function(U, fin, type="log",ylab=FALSE,size=15, plot_type="heatmap")
{
 poly <- getPolygenicity(U) 
 min <- min(c(poly$Me),c(poly$Mt))
 max <- max(c(poly$Me),c(poly$Mt))
 poly %<>% dplyr::filter(Factors %in% fin) %>% mutate(Factors=factor(Factors, levels=rev(fin))) %>%
   pivot_longer(cols=c("Me", "Mt"),values_to = "Poly_score")
 
 if(plot_type=="barplot")
 {
   
   poly$Factors <- factor(poly$Factors, levels=(fin))
   if( type=="log")
   {
     poly$Poly_score = log10(poly$Poly_score)
   }
   if(type == "prop")
   {
     poly$Poly_score = poly$Poly_score/nrow(U)
   }
   p1 = ggplot(poly %>% dplyr::filter(name == "Me"), aes(x=Factors,y=(Poly_score)) ) + geom_bar(stat="identity", fill="#FF5733") + theme_classic(size) + 
     ylab("Effective polygenicity (Me)")
   p2 = ggplot(poly %>% dplyr::filter(name == "Mt"), aes(x=Factors,y=(Poly_score)) ) + geom_bar(stat="identity", fill="#FF5733") + theme_classic(size) + 
     ylab("Total polygenicity (Mt)")
   return(list("me"=p1, "mt"=p2))
 }

 
 
 
 if(type=="log")
 {
   p <- ggplot(poly, aes(x=name,y=Factors,fill=log10(Poly_score))) + geom_tile() + 
     scale_fill_gradient(low="aliceblue", high="red", limits=c(log10(min),log10(max))) + 
     theme_classic(size)  + theme(axis.title.x=element_blank(),legend.position = "bottom",
                           axis.text.x = element_text(angle=45, vjust=1,hjust=1)) + labs(fill="Polygenicity")+
     scale_x_discrete( labels = c("Effective","Total")) 
 }else
 {
   poly$Poly_score <- poly$Poly_score/nrow(U)
   p <- ggplot(poly, aes(x=name,y=Factors,fill=Poly_score)) + geom_tile() + 
     scale_fill_gradient(low="aliceblue", high="red", limits=c(min/nrow(U),max/nrow(U))) + 
     theme_classic(size)  + theme(axis.title.x=element_blank(),legend.position = "bottom",
                           axis.text.x = element_text(angle=45, vjust=1,hjust=1)) + labs(fill="Polygenicity")+
     scale_x_discrete( labels = c("Effective","Total")) 
 }

   if(!ylab)
   {
     p <- dropYLab(p)
   }
 return(p )

}

dropYLab <- function(p)
{
  p  +  theme(axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
}

dropXLab <- function(p)
{
  p  +  theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
}
#tissue.plotting.dat- a table output by the ldsc viz script

nominateNonZeroTissues <- function(tissue.plotting.dat, fdr=0.05)
{
  as.character(sort(factor(unique((tissue.plotting.dat %>% dplyr::filter(factor_tissue_fdr < fdr))$Source),levels=paste0("F",1:100))))
}

plotNiceTissues <- function(ast.factors,tissue.plotting.dat,fdr_pass=0.05, size=15,ylab=FALSE)
{
  #ast.factors <- c("F1","F2", "F6","F8", "F27")
  #View(scaled.V %>% select(Traits, V1, V2, V6,V8, V27))
  focus_trait_tissue <- tissue.plotting.dat %>% dplyr::filter(Source %in% ast.factors) 
  message("dplyr::filtering by factorspecific tissue FDR, note that adjusted for tissue")
  ggplot(focus_trait_tissue %>% dplyr::filter(Coefficient_P_value < 0.01), aes(x=Name,y=-log10(factor_tissue_fdr), fill=color_code)) +
    geom_bar(stat="identity") +facet_wrap(~Source, scales = "free")
  
  factor_order = paste0("F",100:1)
  tissue.order <- unique((focus_trait_tissue %>% arrange(new_category))$p_tissue)
  pass_tissues <- unique((focus_trait_tissue %>% dplyr::filter(Factor_specific_fdr < fdr_pass))$tissue)
  pass_cat_names <- unique((focus_trait_tissue %>% arrange(color_code) %>% dplyr::filter(Factor_specific_fdr < fdr_pass))$new_category)
  pass <- focus_trait_tissue %>% dplyr::filter(tissue %in% pass_tissues) %>% mutate("coef_out" = ifelse(Factor_specific_fdr < fdr_pass, z_score, 0 )) %>%
    mutate("color_code"=ifelse(Factor_specific_fdr < fdr_pass, color_code,"#FFFFFF"))
  p <- ggplot(pass, aes(x= factor(Source, level = factor_order), y =factor(p_tissue, level = tissue.order), 
                   fill = color_code, alpha = coef_out)) + 
    geom_tile(color = "gray")  + scale_fill_identity(guide = "legend", labels =pass_cat_names) + 
    scale_y_discrete(label=function(x) abbreviate(gsub(x,pattern = "_", replacement = " "), minlength = 20)) + 
    ylab("Tissue type") + xlab("Factor Number") +
    guides(fill=guide_legend(title="Tissue Category")) + theme_classic(size) + labs("alpha" = "Z-score") + 
    coord_flip() + theme(axis.text.x = element_blank(), legend.position = "bottom")
  if(!ylab)
  {
    p <- dropYLab(p)
  }
  p
}




plotNiceTissuesBarplot <- function(ast.factors,tissue.plotting.dat,fdr_pass=0.05, size=15,ylab=FALSE)
{
  #ast.factors <- c("F1","F2", "F6","F8", "F27")
  #View(scaled.V %>% select(Traits, V1, V2, V6,V8, V27))
  focus_trait_tissue <- tissue.plotting.dat %>% dplyr::filter(Source %in% ast.factors) 
  message("dplyr::filtering by factorspecific tissue FDR, note that adjusted for tissue")
  factor_order = paste0("F",1:1000)
  tissue.order <- unique((focus_trait_tissue %>% arrange(new_category))$p_tissue)
  pass_cat_names <- unique((focus_trait_tissue %>% arrange(color_code) %>% dplyr::filter(Factor_specific_fdr < fdr_pass))$new_category)
  pass <- focus_trait_tissue %>% dplyr::filter(Factor_specific_fdr < fdr_pass) %>% rename("Factor"=Source)
  
  count.by.cat <- pass %>% group_by(Factor, category,color_code) %>% summarize("Count"=n()) 
  #append on any missing factors
  count.by.cat <- rbind(count.by.cat, 
                        data.frame("Factor"=ast.factors[!(ast.factors %in% count.by.cat$Factor)], category=NA, color_code="#FFFFFF", Count=0)) %>%
    mutate("Factor"=factor(Factor, levels = factor_order )) 
  
         p <- ggplot(count.by.cat, aes(x=Factor, y= Count, fill=color_code)) + geom_bar(stat="identity") + 
           scale_fill_identity(guide = "legend", labels =pass_cat_names) + theme_classic(size) + theme(legend.position = "bottom")+
           labs(fill="")+ ylab("# tissue marker enrichments")

         if(!ylab)
         {
           p <- dropYLab(p)
         }
         p
}







## Calculate scores based on p-values
#getPleiotropyScores(p.vals.all,full.U, threshold=1e-5)

getPleiotropyScores <- function(pvals,U.mat, threshold=1e-5)
{
  
  pleio.scores <- data.frame("SNPs" =pvals[,1], "pleio" = apply(pvals[,-1],1,function(x) sum(x < threshold, na.rm=TRUE)/sum(!is.na(x)))) %>%
    set_colnames(c("SNPs", "pleio")) %>% dplyr::filter(SNPs %in% U.mat$SNPs)

  left_join(U.mat, pleio.scores, by="SNPs")

}

#Full U has SNPs in teh front and then all the factors
getTopAndBottomSNPSByFactor <- function(full.U,u.with.plieo)
{
  top.snps.by.factor <- lapply(2:ncol(full.U), function(i) which(full.U[,i]^2 > quantile(full.U[,i]^2,probs = 0.99)))
  bottom.snps.by.factor <- lapply(2:ncol(full.U), function(i) which(full.U[,i]^2 < quantile(full.U[,i]^2,probs = 0.995)))
  
  names(top.snps.by.factor) <- colnames(full.U)[-1];names(bottom.snps.by.factor) <- colnames(full.U)[-1]
  #Top snps
  scores.by.factor <- lapply(top.snps.by.factor, function(x) u.with.plieo$pleio[x])
  pseudo.pleio.df <- do.call("rbind", lapply(names(scores.by.factor), function(n) data.frame(scores.by.factor[[n]], "source"=n))) %>% 
    set_colnames(c("Pleioscore", "factor")) %>% mutate("where"="Top 1%")
  
  scores.by.factor.lower <- lapply(bottom.snps.by.factor, function(x) u.with.plieo$pleio[x])
  pseudo.pleio.bottom.df <- do.call("rbind", lapply(names(scores.by.factor.lower), function(n) data.frame(scores.by.factor.lower[[n]], "source"=n))) %>% set_colnames(c("Pleioscore", "factor"))%>% mutate("where"="Bottom")
  rbind(pseudo.pleio.bottom.df,pseudo.pleio.df)
}

#Factors of interest- then ames of teh factors to pick out (U)
#joined.pleio- a list of the pleiotropy score of each SNP, the SNP's factor, and where it falls (either top1% or bottom.)
sidewaysPleioBarplot <- function(joined.pleio,factors_of_interest,size=15, ylab=FALSE)
{
  
  #ggplot(joined.pleio %>% dplyr::filter(factor %in% factors_of_interest), aes(x=factor,y=Pleioscore/137, fill=where)) + geom_boxplot() + theme_classic()
  #joined.pleio %>% dplyr::filter(factor %in% factors_of_interest) %>% group_by(factor,where) %>% summarize("count"=n(), "avg"=mean(Pleioscore))
  
  
  #ggplot(joined.pleio %>% dplyr::filter(factor %in% factors_of_interest), aes(x=factor,y=Pleioscore, fill=where)) + geom_violin() + theme_classic()
  joined.pleio %<>% dplyr::filter(factor %in% factors_of_interest) %<>% mutate("factor"=factor(factor,levels=rev(factors_of_interest)))
  p <- ggplot(joined.pleio, aes(x=factor,y=Pleioscore, fill=where)) + geom_boxplot() + theme_classic(size) + coord_flip() +
    theme(legend.position="bottom") + guides(fill=guide_legend(ncol=1)) + 
    labs(fill=bquote(U^2)) + xlab("Pleiotropy score")
  if(!ylab) p <- dropYLab(p)
  p
}


plotSelectivePressure <- function(f.in, 
                                  select.dat="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/selective_pressure/s_scores.tsv",
                                  ylab=FALSE,size=15, type = "heatmap")
{
  selection <- fread(select.dat) 
  min.s <- min(selection$S_hat)
  max.s <- max(selection$S_hat)
  selection %<>% dplyr::filter(Factor %in% f.in)
  selection$Factor <- factor(selection$Factor, levels=rev(f.in))
  if(type=="heatmap")
  {
    p <- ggplot(selection, aes(y=Factor, fill=S_hat,x=1)) + geom_tile() +
      scale_fill_gradient(high="bisque",low="blue",limits = c(min.s,max.s)) + theme_classic(size) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),legend.position = "bottom") + labs(fill=bquote(hat(S)))
  }
  if(type=="barplot")
  {
    selection$Factor <- factor(selection$Factor, levels=(f.in))
    p <- ggplot(selection,(aes(x = Factor, y=S_hat))) + geom_bar(stat="identity", fill="#33AFFF") + theme_classic() + 
      ylab("Selection coef")
  }

  if(!ylab) p <- dropYLab(p)
  p
  

}

#sept 7 dropped pleioitropy
#overallNicePlot <- function(eval.factors,trait, scaled.V,U,joined.pleio,
#                            tissue.dat = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/NONE_ldsc_enrichment_Multi_tissue_chromatin/plotting_table.csv",
#                            ...)
overallNicePlot <- function(eval.factors,trait, scaled.V,U,select_dat,
                           tissue.dat = "/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_61K/NONE_ldsc_enrichment_Multi_tissue_chromatin/plotting_table.csv",
                           ...)
{
  tissue.plotting.dat <- fread(tissue.dat)
  barplots <- cowplot::plot_grid(plotlist = traitBarplots(trait, scaled.V,factor_list=eval.factors),ncol = 1)
  tissues <- plotNiceTissues(paste0("F",eval.factors), tissue.plotting.dat,...)
   #############
  #pleio.plot <- sidewaysPleioBarplot(joined.pleio,paste0("U",eval.factors))
  
  ############ Selective pressure
  selective.pressure <- plotSelectivePressure(paste0("U",eval.factors),select.dat = select_dat)
  
  ########### Number of effective SNPs and  sample
  polygen <- plotPolygenicity(U, paste0("U",eval.factors),type = "prop")
  #return(list(barplots,tissues, pleio.plot,selective.pressure,polygen))
  return(list(barplots,tissues, selective.pressure,polygen))
}


######
agreementWithLDSC <- function(V,traits, path="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/panUKBB_complete/rg_ldsr/tabular/")
{
  joined.all <- do.call("rbind",lapply(list.files(path, pattern = "*rg_report.tsv"), function(x) fread(paste0(path,x))))
  joined.all$trait1 <- stringr::str_extract(string = joined.all$p1, pattern="EUR_ONLY\\.(.*)\\.sumstats\\.gz",group = 1) %>% 
    gsub(x=.,pattern="non-albumin", replacement="non.albumin")%>% gsub(x=.,pattern="^_", replacement="X_")
  joined.all$trait2 <- stringr::str_extract(string = joined.all$p2, pattern="EUR_ONLY\\.(.*)\\.sumstats\\.gz",group = 1) %>% 
    gsub(x=.,pattern="non-albumin", replacement="non.albumin") %>% gsub(x=.,pattern="^_", replacement="X_")
  paired.entries <- lapply(apply(V,2,function(x) which(x!=0)), function(x) traits[x])
  names(paired.entries) <- paste0("F",1:ncol(V))
  all.ps <- c()
  all.rgs <- NULL
  all.dat <- NULL
  factor.dat <- list()
  for(i in 1:length(paired.entries))
  {
    factor <- paste0("F",i)
    l <- unlist(paired.entries[[i]])
    
    if(length(l) == 1)
    {
      next
    }
    all.in <- dplyr::filter(joined.all, trait1 %in% l & trait2 %in% l)
    if(nrow(all.in) < 1)
    {
      print(i)
      message("BAd man")
      print(l)
      next
    }
    #all.ps <- c(all.ps, all.in$p)
    #all.rgs <- c(all.rgs, all.in$rg)
    all.dat <- rbind(all.dat, data.frame("Factor"=factor,"rg"=all.in$r, "p"=all.in$p, "n"=nrow(all.in)))
    #factor.dat[[factor]] <- list("rg"=all.in$rg,"p"=all.in$p, "n"=nrow(all.in) )
    #I want the number of pairs, the distribution on Rg, and the distribution on P
  }
  all.dat$Factor <- factor(all.dat$Factor, levels = paste0("F",1:100))
  #This will always be biased towards HM3, because that is what we consider over. derp.
  f <- ggplot(all.dat, aes(x=Factor,y=-log10(p+1e-20))) + geom_boxplot() + theme_bw() + geom_hline(color="blue", yintercept = -log10(0.01)) + 
    xlab("Factor") + ylab("-log10(p) vals") + ggtitle("Distribution of p-values per factor")
  s <- ggplot(all.dat, aes(x=Factor,y=rg^2)) + geom_boxplot() + theme_bw() + xlab("Factor") + ylab(bquote(r[g]^2)) + ggtitle("Distribution of rg^2 per factor")
  return(list(all.dat, f,s))
  
}


##Pleioscores


### Scores on blieotropy
#ret$U, u.with.plieo$pleio, maf=u.with.plieo$maf)
pleioScoreByWilcox <- function(input.factor, pleio_scores, top_perc = 0.99, alpha=0.01, ntests=100, maf=NA, pleio_dir="greater", diff_dat=TRUE,
                               method="wilcox")
{
  # Find indices of top and bottom SNPs based on the top_perc
  top.snps.by.factor <- which(input.factor^2 >= quantile(input.factor^2,probs = top_perc))
  bottom.snps.by.factor <- which(input.factor^2 < quantile(input.factor^2,probs = top_perc))
  bottom.snps.by.factor_log <- input.factor^2 < quantile(input.factor^2,probs = top_perc)
  pleio.top <- pleio_scores[top.snps.by.factor]
  count.below.thresh <- 0
  ret_vals <- c()
  #build group thresholds if within maf window:
  if(any(!is.na(maf)))
  {
    maf.cuts <- cut(maf, breaks=10)
    maf.upper <- table(maf.cuts[top.snps.by.factor])
    maf.lower <- table(maf.cuts[bottom.snps.by.factor])
    if(any(maf.lower < maf.upper))
    {
      message("Warning- fewer in bg group than in selected group. Be careful with results.")
    }
  }
  
  for(n in 1:ntests)
  {
    if(any(!is.na(maf)))
    {
      pleio.rand <- c()
      for(i in 1:length(maf.upper))
      {
        l=names(maf.upper)[i]
        l_val=maf.upper[i]
        #get pleio scores in the same bin
        pleio.rand <- c(pleio.rand, sample(x=pleio_scores[(maf.cuts == l) & (bottom.snps.by.factor_log)],
               size=l_val))
        
      }
    }else
    {
      random.equal.length <- sample(x=bottom.snps.by.factor, size=length(top.snps.by.factor))
      pleio.rand <- pleio_scores[random.equal.length] #since going by index, this works.
    }
    
    #Now, to perform the test. There are many options:
    #diff_dat specieis that you actually want the statistics themsleves returned, and not just if they are at or above some threhsold.
    if(!diff_dat)
    {
      get.stat <- wilcox.test(pleio.top,pleio.rand,alternative=pleio_dir)
    }else if(diff_dat & method == "t.test") #Do a t-test, get those stats
    {
      get.stat_all <- t.test(pleio.top,pleio.rand,alternative=pleio_dir)
      get.stat<- get.stat_all$statistic
    }else if(method=="fold_change" & diff_dat) #currently underdeveloped- just look at the ratio. Not obvious how to get SE bars on this.
    {
      get.stat <- mean(pleio.top)/mean(pleio.rand)
    }
    else if(method=="logit" & diff_dat) #currently underdeveloped- just look at the ratio. Not obvious how to get SE bars on this.
    {
      outcome = c(rep(1,length(pleio.top)),rep(0,length(pleio.rand)))
      pleio.scores <- c(pleio.top,pleio.rand)*137
      logit_vs <- glm(outcome ~ pleio.scores, family = "binomial")
      get.stat <- exp(coef(logit_vs))[2]
    }
    else
    {
      get.stat <- wilcox.test(pleio.top,pleio.rand,alternative=pleio_dir, conf.int = TRUE)$estimate
    }
    
    if(method == "t.test" | method == "logit"| method== "fold_change" )
    {
      ret_vals <- c(ret_vals,get.stat)
    }
    else if( (!diff_dat) & (get.stat$p.value < alpha)  )
    {
      count.below.thresh <- count.below.thresh + 1
    }else
    {
      ret_vals <- c(ret_vals,get.stat)
    }
  }
  if((!diff_dat))
  {
    return(count.below.thresh)
  }else
  {
    return(ret_vals)
  }
  
}

loadGeneEnrichments <- function(factor, type, 
                                path="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_final/top_elements_by_factor/")
{
  ret.list <- list()
  old.mapping <- factor
  dir=paste0(path,type)
  ret.list$all_genes <- fread(paste0(dir, "/F", old.mapping, "_snp_gene_rank.csv"))
  ret.list$set_genes <- fread(paste0(dir, "/F", old.mapping, "_genes_in_factor.txt"),header = FALSE)$V1
  ret.list$bg_genes <- fread(paste0(dir, "/F", old.mapping, "_background_genes.txt"),header = FALSE)$V1
  
  ret.list$enrichments <- fread(paste0(dir, "/gene_factor_enrichments/F", old.mapping, ".joined_enrichments.txt"))
  ret.list$enrichments_unique <- fread(paste0(dir, "/gene_factor_enrichments_unique/F", old.mapping, ".joined_enrichments.txt"))
  ret.list
  
}

pleioScoresAll <- function(u.dat, pleio_scores,...)
{
  apply(u.dat,2,function(u) pleioScoreByWilcox(u,pleio_scores,...))
}
####
plotPVEPretty <- function(pve_vector, factors, cin="gray", size=15)
{
  to.plot <- data.frame(factors,pve_vector) %>% mutate(factors = factor(factors, levels=factors))
  ggplot(to.plot, aes(x=factors,y=pve_vector))+geom_bar(stat="identity", color=cin) + theme_classic(size) + ylab("PVE") + 
    xlab("Factor")
}



## eenrichment test:
#testEnrichment(f11.tp$set_genes, factor_titles, c(late_mk,mk_hmt),conf.level=0.90,alternative="greater")
testEnrichment <- function(set_genes, bg_genes, gene_set,conf.level=0.95,alternative="two.sided")
{
  contingency.table <- t(matrix(c(sum(set_genes %in% gene_set),sum(!(set_genes %in% gene_set)),
                                  sum(bg_genes %in% gene_set),sum(!(bg_genes %in% gene_set))), nrow = 2))
  list("test" = fisher.test(contingency.table,alternative=alternative,conf.int =conf.level ), "tab"=contingency.table)
}

#testEnrichment(unique(f4.tp$set_genes),unique(f4.tp$bg_genes), early_Megakaryopoiesis),
testEnrichmentPermutedUnique <- function(set_genes, bg_genes, gene_set,n_permute=10000, seed=NA,...)
{
  if(!is.na(seed)){set.seed(seed)}
  permute.list <- list()
  n.samp <- length(set_genes)
  for(i in 1:n_permute)
  {
    #svMisc::progress(i,max.value=n_permute)
    subsample.genes <- sample(bg_genes, n.samp)
    permute.list[[i]] <- testEnrichment(unique(subsample.genes), unique(bg_genes), gene_set,...)
    if(i %% 10000 == 0)
    {
      message("Iteration ", i)
    }
  }
  #Now on the real set
  real.vals <- testEnrichment(unique(set_genes), unique(bg_genes), gene_set,...)
  #Interpret this elsewhere.....
  perm.ors <- sapply(permute.list, function(x) x$test$estimate)
  list("set_score"=real.vals, "permuted_p"=sum(perm.ors >= real.vals$test$estimate)/n_permute, "n_permute"=n_permute, "permuted_data"=permute.list)
}


#factor_sets = list(f4.tp$set_genes, f11.tp$set_genes, f23.tp$set_genes)
#factor_titles <- c("F4","F11","F23")
#sets_to_test <- list(early_Megakaryopoiesis,late_and_proto,cytoskeleton,platelet_function)
#factor_sets,full.bg,factor_titles, sets_to_test,n_permute=1000
tabEnrichmentTestPermuted <- function(factor_sets,full.bg,factor_titles, sets_to_test,gene_set_names,n_permute=1000,...)
{
  all.tests <- list() 
  for(i in 1:length(factor_sets))
  {
    message("Factor set ",i, " of ", length(factor_sets) )
    f <- factor_sets[[i]]
    fname <- factor_titles[[i]]
    all.tests[[factor_titles[[i]]]] <- lapply(sets_to_test, function(set)  testEnrichmentPermutedUnique(f,full.bg, set, n_permute = n_permute,...))
  }

  #Build the table
  do.call("rbind", lapply(1:length(all.tests), function(i) {
    dat <- all.tests[[i]]; fname = factor_titles[i]
    data.frame("Factor" = rep(fname,length(dat)), "p-val" = sapply(dat, function(x) x$permuted_p), 
               "OR"=  sapply(dat, function(x) x$set_score$test$estimate), 
               "LOWER"=sapply(dat, function(x) x$set_score$test$conf.int[1]),
               "UPPER"=sapply(dat, function(x) x$set_score$test$conf.int[2]),"Phase"=gene_set_names,
               "num_genes"=sapply(dat, function(x) x$set_score$tab[1]),
               "set_size"=sapply(sets_to_test, length),
               "true_test_p"=sapply(dat, function(x) x$set_score$test$p.value))
  })) %>%
    mutate("p.val" = ifelse(p.val == 0, 1/n_permute,p.val)) %>% mutate("bonf.p" = p.adjust(p.val, method="bonferroni")) %>%
    mutate("FDR" = p.adjust(p.val, method="BH")) %>%
    mutate("p_sig"=ifelse(bonf.p < 0.05, "p < 0.05", "p > 0.05"))
}


loadMAFData <- function()
{
  full.snps <- fread("/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz")
  full.snps <- filter(full.snps,rsid %in% ret$snp.ids) #There are 2 redundant SNPs it looks like.
  full.snps$maf <- sapply(full.snps$af_EUR, function(x) ifelse(x > 0.5, 1-x, x))
  full.snps
}


evalTissuesLikeEQTLs <- function(tiss.df, factor_id)
{
  tiss.df %>%  filter(Source == factor_id) %>% group_by(new_category) %>% mutate("category_bfp"=p.adjust(Coefficient_P_value,method = "bonferroni")) %>% slice_min(category_bfp) %>%
    ungroup() %>% mutate("max_tissue_fdr_adjusted"=p.adjust(category_bfp, method="BH"))

  
}

#New functions for new simulations:
#source("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/R/data_decorrelation.R")
library(gleanr)

#Functions for decomposing string identifiers. Very bespoke.
getScore <- function(strin)
{
  sum(gsub(x=strin, pattern="overlap-",replacement = "") %>% stringr::str_split(., pattern ="-", n=3) %>%
        sapply(., as.numeric))
}
numWithOverlap <- function(strin)
{
  sum(gsub(x=strin, pattern="overlap-",replacement = "") %>% stringr::str_split(., pattern ="-", n=3) %>%
        sapply(., as.numeric) != 0)
}


#VISUALIZE THE EIGENVALUE PLOTS..
eigenAnalysis <- function(pca.sims, file.list.betas)
{
  
  eigen.df <- NULL
  for(i in 1:length(pca.sims))
  {
    eigen.df <- rbind(eigen.df, data.frame("d"=pca.sims[[i]]$d, 
                                           "case" = gsub(file.list.betas[[i]],pattern = ".BETA.csv", replacement = ""), 
                                           "eig_rank"=1:length(pca.sims[[i]]$d),
                                           "pve"=pca.sims[[i]]$d^2/sum(pca.sims[[i]]$d^2)))
  }
  eigen.df <- eigen.df %>% separate(case, sep = "\\.", into = c("seed", "overlap", "replicate")) %>% rowwise() %>%
    mutate("overlap_score"=getScore(overlap))%>% ungroup() %>% arrange(overlap_score)
  
  two <- ggplot(eigen.df, aes(x = as.factor(overlap_score), y=as.numeric(d))) + geom_boxplot(varwidth = TRUE) + xlab("Sum of overlapping samples") + ylab("Eigenvalues") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(list(eigen.df, two))
}

tableSimResults <- function(len, gd.point, kaiser.point, elbow.point, flashr.runs, gleaner.runs.std, gleaner.runs.covar, files.betas)
{
  factors.by.method <- data.frame("gd"=unlist(gd.point)[1:len],"kaiser"=kaiser.point[1:len], "elbow"=unlist(elbow.point)[1:len], 
                                    "flash"=sapply(flashr.runs, function(x)x$nfactors),
                                    "gleaner_noadj" = unlist(sapply(gleaner.runs.std, function(x) ifelse(ncol(as.matrix(x$V)) & all(x$V==0),0,1))),
                                    "gleaner_covar" = unlist(sapply(gleaner.runs.covar, function(x) ifelse(ncol(as.matrix(x$V)) & all(x$V==0),0,1))),
                                    "case" = sapply(files.betas[1:len], function(x) gsub(x,pattern = ".BETA.csv", replacement = ""))) %>% 
    separate(case, sep = "\\.", into = c("seed", "overlap", "replicate")) %>% rowwise() %>%
    mutate("overlap_score"=getScore(overlap), "study_with_overlap"=numWithOverlap(overlap)*3) %>%
    ungroup() 
  
  
  factors.by.method.long <- factors.by.method %>% pivot_longer(cols=c("gd", "kaiser", "elbow", "flash", "gleaner_noadj", "gleaner_covar"),names_to = "method")
  return(list("std"=factors.by.method, "long"=factors.by.method.long))
  
}

tableSimResultsFromList <- function(l)
{
  possible.columns.to.pars <- c("gleaner.runs.covar","gleaner.runs.std",
                                  "whitened_sim.pca.std","whitened_sim.pca.alt","whitened_gd.point","whitened_elbow.point","whitened_kaiser.point",
                                  "flashr.runs.se","elbow.point","gd.point","sim.pca.alt","sim.pca.std", "kaiser.point")
  df.tb <- list()
  len=length(l$gleaner.runs.std)
  for(p in possible.columns.to.pars)
  {
    if(p %in% names(l))
    {
      if(p %in% c("gd.point", "kaiser.point", "elbow.point", "whitened_gd.point","whitened_elbow.point","whitened_kaiser.point"))
      {
        df.tb[[p]] <- unlist(l[[p]])[1:len]
      }
      if(p %in% c("gleaner.runs.std", "gleaner.runs.covar"))
      {
        df.tb[[p]] <- unlist(sapply(l[[p]], function(x) ifelse(all(x$V==0),0,ncol(as.matrix(x$V)))))
      }
      if(p == "flashr.runs.se")
      {
        df.tb[[p]] <- sapply(l[[p]], function(x)x$n_factors)
      }
    }
  }
  df.tb[["case"]] = sapply(l$files.betas[1:len], function(x) gsub(x,pattern = ".BETA.csv", replacement = ""))
  
  df <- data.frame(do.call("cbind", lapply(df.tb, function(x) x))) %>% 
    separate(case, sep = "\\.", into = c("seed", "overlap", "replicate")) 
  df$overlap_score = sapply(df$overlap, function(o) getScore(o))
  df$study_with_overlap = sapply(df$overlap, function(o) numWithOverlap(o)*3)
    #mutate("overlap_score"=getScore(overlap), "study_with_overlap"=numWithOverlap(overlap)*3) %>%
    #ungroup() 
  pivot.cols <-colnames(df)[1:9]
  factors.by.method.long <- df %>% pivot_longer(cols=all_of(pivot.cols),names_to = "method")
  return(list("std"=df, "long"=factors.by.method.long))
  
}



##Whitening tools
estimate_null_correlation_simple <- function(z, z_thresh = 2, est_cor = TRUE)
{
  #assumes z input
  max_absz = apply(abs(z), 1, max)
  nullish = which(max_absz < z_thresh)
  if (length(nullish) < ncol(z)) {
    stop("not enough null data to estimate null correlation")
  }
  nullish_z = z[nullish, ]
  Vhat = cov(nullish_z)
  if (est_cor) {
    Vhat = cov2cor(Vhat)
  }
  return(Vhat)
}

#the whole thing.
#z is the raw data data.
blockCorrelationMatrix <- function(z,cor_thresh = 0.2, blocks = NULL)
{
  est <- estimate_null_correlation_simple(z)
  if(any(is.na(blocks)))
  {
    blocks <- create_blocks(est, cor_thr=cor_thresh)
  }

  ret <- blockify(est, blocks)
  if(!isSymmetric(ret))
  {
    print("not symmetric?")
    print(ret)
  }
  return(ret)
}

#Functions to help process
#Get all GWAS
#Here doing as z-scores
readInGWAS <- function(sim.results, dir)
{
  gwas.list <- list()
  for(i in 1:length(sim.results))
  {
    x <- sim.results[i]
    gwas.list[[i]] <- fread(paste0(dir,x))
  }
  return(gwas.list)
}

#Get PCA for all
pcaOnAll <- function(gwas.list, cols = "ALL")
{
  if(cols == "ALL")
  {
    lapply(gwas.list, function(x) svd(scale(x)))
  }else
  {
    lapply(gwas.list, function(x) svd(scale(x[,cols]))) #limit to just certain ones that contiain the data of interest.
  }

}

GLEANEROnAll <- function(gwas.list, cols = "ALL", covar=TRUE)
{
  setwd("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/"); load_all()
  if(covar)
  {
    #compare runs to choose the best lambda
    X <- gwas.list
    message("Currently just using Z-scores, need to update this...")
    W <- matrix(1, nrow = nrow(X), ncol = ncol(X))
    C <- truncatedCovarMatrix(X, z_thresh = 2)
    
  }else
  {
    lambda = 1
    sapply(paste0("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gleanr/R/", files.source), source)
    res <- gleaner(X,W, paste0("rs", 1:nrow(X)), paste0("T",1:ncol(X)), K=K,C=C,...)
  }
}



#Assess spearman correlation for all
  #spearmanVsPCs(true.null.pcs, c(0,0,0,1,1,1),trait.names = c("A1", "B1","C1", "A2", "B2", "C2") )
## spearmanVsPCs(true.null.sims,true.null.pcs, c(0,0,0,1,1,1),trait.names = c("A1", "B1","C1", "A2", "B2", "C2"),stat.m = "kendall")

spearmanVsPCs <- function(run.name, pc.list, spearman.key,trait.names = c("A1", "B1", "A2", "B2"), stat.m = "spearman", nsv="kaiser")
{
  rho=list()
  p=list()
  for(i in 1:length(pc.list))
  {
    svs <- pc.list[[i]]
    pve <- svs$d^2/sum(svs$d^2)
    rownames(svs$v) <- trait.names
    #want to update this- looking at average
    if(nsv!="kaiser")
    {
      sig.pves <- 1:as.numeric(nsv)
    } else
    {
      sig.pves <- which(pve > mean(pve)) #Kaiser based method.
    }
  
    #sig.pves <- which(pve > (1/length(svs$d)))
    rho[[i]] <- lapply(sig.pves, function(x) cor.test(y=spearman.key, x = svs$v[,x], method = stat.m)$estimate)
    p[[i]] <- lapply(sig.pves, function(x) cor.test(y=spearman.key, x = svs$v[,x], method = stat.m)$p.value)
  }
  ret <- NULL
  for(j in 1:length(rho))
  {
    for(i in 1:length(rho[[j]]))
    {
      ret <- rbind(ret, c(gsub(run.name[[j]],pattern = ".gwas.csv", replacement = ""), i, as.numeric(rho[[j]][[i]]), as.numeric(p[[j]][[i]])))
    }
  }
    colnames(ret) <- c("run", "PC", "rho", "pval")
    ret.frame <- data.frame(ret) %>% separate(run, sep = "\\.", into = c("seed", "rep", "overlap"))
    #ret$overlap <- as.numeric(ret$overlap);ret$rho <- as.numeric(ret$rho);ret$pval <- as.numeric(ret$pval)
  return(ret.frame)
}

# Reconstruction error




# Pve change plot
pveChangePlot <- function(pc.list, list.names)
{
  null.pve <- do.call("rbind", lapply(1:length(pc.list), function(i) data.frame("pve" = pc.list[[i]]$d^2/sum(pc.list[[i]]$d^2)) %>% mutate("sim_id" = gsub(pattern=".gwas.csv", replacement = "",x = list.names[i])) %>%
                                        tidyr::separate(sim_id, into=c("seed","sim_num", "overlap"), sep = "\\.", remove = TRUE) %>% mutate("pc"=1:length(pc.list[[i]]$d)))) %>%
    mutate("perc_overlap" = round(as.numeric(overlap)/15000, digits = 2))
  ret.plot <- ggplot(null.pve, aes(x = as.factor(perc_overlap), y=pve )) + geom_boxplot() +
    geom_smooth(method = "lm", se=FALSE, color="blue", aes(group="1")) +
    facet_wrap(~pc,scales = "free_y" ) + theme_classic(15) + xlab("Percent overlap") + ylab("PVE")
  return(list("df" = null.pve, "plot" = ret.plot))
}

#Pairwise R2 plot:
#vp.null.sims
#look at 1-6
#lapply(1:6, function(x) cor(vp.null.zscores[[x]])[c])

pairwiseR2DF <- function(z.scores,sim.names, combs = matrix(c(c("A1","A1", "B1", "A2", "A2", "B2"),c("B1","C1", "C1", "B2", "C2", "C2")), ncol = 2))
{
  all.cors <- lapply(z.scores, function(x) cor(x)^2)
  combs.vect <- apply(combs, 1, function(x) paste0(x[1], ":", x[2]))
  tab.combs <- lapply(all.cors, function(c)  cbind(apply(combs, 1, function(x) paste0(x[1], ":", x[2])), c[combs]))
  joined.combs <- do.call("rbind", lapply(1:length(tab.combs), function(i) data.frame(tab.combs[[i]], "rep" = gsub(sim.names[i],pattern = ".gwas.csv", replacement = "")))) %>%
    set_colnames(c("GWAS_entries", "R2", "run")) %>%
    separate(run, sep = "\\.", into = c("seed", "rep", "overlap")) %>%
    mutate("cohort" = substr(GWAS_entries,start = 2, stop = 2) )
  joined.combs
}


#KS test tools
condKSTest <- function(dat, fcol, scol, thresh = 0.1)
{
  cond.tests <- dat[,scol][dat[,fcol] < thresh]
  ks.test(cond.tests, y="punif")
}

multiCondKS <- function(dat, pairs, pval = TRUE)
{
  if(pval) {return(lapply(pairs, function(x) condKSTest(dat, x[1], x[2])$p.value))}
  else
  {
    return(lapply(pairs, function(x) condKSTest(dat, x[1], x[2])$statistic))
  }
}


ksTestTab <- function(zin, sim.names, pairs = list(c("A1", "C1"), c("B1", "C1"), c("A2", "C2"), c("B2", "C2")), pair.names=c("A1:C1", "B1:C1", "A2:C2", "B2:C2"), pval = TRUE)
{
    ks.dat <- NULL
  for(tab in zin)
  {
    #make it p-values
    pvals <- data.frame(apply(abs(tab), 2, function(x) pnorm(-x)*2))
    if(pval)
    {
      std.unif <- apply(pvals, 2, function(x) (ks.test(x, y = "punif"))$p.value)
    } else{
      std.unif <- apply(pvals, 2, function(x) (ks.test(x, y = "punif"))$statistic)
    }

    conditionals <- unlist(multiCondKS(pvals, pairs, pval = pval))
    ks.dat <- rbind(ks.dat, c(std.unif, conditionals))
  }

  ks.df <- data.frame(ks.dat) %>% set_colnames(c(colnames(tab), pair.names)) %>%
    mutate("run" = gsub(sim.names,pattern = ".gwas.csv",replacement = "")) %>% pivot_longer(cols = all_of(pair.names), "conditional_test") %>% separate(run, sep = "\\.", into = c("seed", "rep", "overlap"),remove = FALSE)
  ks.df$overlap <- as.numeric(ks.df$overlap)
  ks.df
}


#R2 performacne testing...
#pcas.each <- signal.simple.pcs
#sim.names <- signal.simple.sims
testFactorizationPerformance <- function(sim.names,pcas.each,npcs="detect", whitening = c())
{
  baseline <- which(grepl(sim.names, pattern = "\\.0.gwas.csv"))
  base.names <- sim.names[baseline]
  perf.tab <- NULL
  for(s in sim.names)
  {
    if(s %in% base.names)
    {
      message("Baseline run, skippping")
    }else
    {
      query <- str_split(s, pattern = "\\.")[[1]]
      base.match <- paste0(query[[1]],".", query[[2]], ".0.gwas.csv") #which is the correct one we want?
      base.pca <- pcas.each[[which(sim.names == base.match)]]
      query.pca <- pcas.each[[which(sim.names == s)]]
      if(npcs == "detect")
      {
        max.pc <- max(which(base.pca$d > mean(base.pca$d)))
      }else
      {
        max.pc = npcs
      }
      max.r2 <- evaluteFactorConstruction(as.matrix(base.pca$u[,1:max.pc]), as.matrix(base.pca$v[,1:max.pc]), as.matrix(query.pca$u[,1:max.pc]), as.matrix(query.pca$v[,1:max.pc]))
      if(length(whitening) > 0)
      {
        query.wht <- whitened.pcas.each[[which(sim.names == s)]]
        query.wht.raw <- whitened.pcas.raw[[which(sim.names == s)]]
        max.r2.wht <- evaluteFactorConstruction(as.matrix(base.pca$u[,1:max.pc]), as.matrix(base.pca$v[,1:max.pc]), as.matrix(query.wht$u[,1:max.pc]), as.matrix(query.wht$v[,1:max.pc]))
        max.r2.wht.raw <- evaluteFactorConstruction(as.matrix(base.pca$u[,1:max.pc]), as.matrix(base.pca$v[,1:max.pc]), as.matrix(query.wht.raw$u[,1:max.pc]), as.matrix(query.wht.raw$v[,1:max.pc]))
        perf.tab <- rbind(perf.tab, c(query[[1]], query[[2]], query[[3]], max.r2.wht, "pca_wht"))
        perf.tab <- rbind(perf.tab, c(query[[1]], query[[2]], query[[3]], max.r2.wht.raw, "pca_wht_raw"))
      }
      perf.tab <- rbind(perf.tab, c(query[[1]], query[[2]], query[[3]], max.r2, "pca_std"))

    }
  }
  data.frame(perf.tab) %>% set_colnames(c("seed", "rep", "overlap", "L_R2", "F_R2", "method")) %>% mutate("overlap" = as.numeric(overlap))
}


#forthe factor based ones:
mergeTabularResults <- function(query, base.dir = "/scratch16/abattle4/ashton/snp_networks/scratch/cohort_overlap_exploration/simulating_factors/custom_easy/simulation_outputs/", scale = TRUE)
{
  library(magrittr)
  library(dplyr)
  library(data.table)
  library(stringr)
  d <- list.dirs(base.dir,recursive = FALSE)
  t <- NULL
  library(tidyr)
  file.ext = "/factorization_results//summary.tabular_performance.tsv"
  if(!scale){
    file.ext="/factorization_results//summary.noscale.tabular_performance.tsv"
  }
  for(p in d)
  {
    if(grepl(query, p))
    {
      maf <- stringr::str_extract(string =p,pattern = "MAF-([\\w\\d]+)_eur",group=1)
      if(file.exists(paste0(p, file.ext)))
      {
        print(p)
        n.samp <- stringr::str_extract(string =p,pattern = "N-([0-9e\\.\\+\\w]+)_RHO", group=1)
        c <- fread(paste0(p, file.ext)) %>% 
          mutate("MAF" = gsub(pattern = "maf",replacement = "",x=maf), "N"= n.samp) %>% mutate("fname"=basename(p))
        if(length(c) > 0)
        {
          t <- rbind(t, c)
        }
      }else
      {
        message("Incomplete or terminated run- tabular performance summary missing.")
        message(p)
      }
      
    }
  }
  t #%>% mutate("MAF" = as.numeric(MAF))# %>% mutate("noise" = 1/(MAF*2*(1-MAF) * as.numeric(N)))
}

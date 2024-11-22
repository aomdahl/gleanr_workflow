#Generate phenotypes with specific relationships and patterns and noise, based on input genotype information and SNP weights
#Returns a list of phenotypes of length nrow(genotype) with a little noise added on

#snp.weights- which snps are active?
#build your phenotypes using specified snps and genotypes
#phenotypeBuilder(snp.weights, pop1.geno[,1:n.snps], 1, type = "1+2=3")
phenotypeBuilder <- function(snp.weights, genotype, env_component, type = "1+2=3")
{
  message("genotypes are unscaled- recommend scaling them to have mean 0, variance 1, so we can get a better grip on heritability")
  message("Input environmental terms: ", paste(env_component))
  print(dim(genotype))
  print(length(snp.weights[[1]]))
  ret.pheno = list()
  if(length(env_component) ==1)
  {
	  env_component <- rep(env_component,3)

  }

  if(is.null(snp.weights) & type == "1+2=3")
  { 
    message("Warning- noise on phenotypes may be more than you think....")
    ret.pheno[[1]] =  rnorm(nrow(genotype),2,1) 
    ret.pheno[[2]] =  rnorm(nrow(genotype),2,1)
    ret.pheno[[3]] = 0.5*ret.pheno[[1]] +  0.5*ret.pheno[[2]]
  }else if(is.null(snp.weights) & type == "bmi_like"){
    #Numbers selected from UKBB-based parameters; see script ______
    ret.pheno[[1]] = rnorm(nrow(genotype), mean= 169, sd= 9.3)/100 #scaled from cm to m
    ret.pheno[[2]] = -80 + 0.93 *(ret.pheno[[1]]*100)  + rnorm(nrow(genotype),sd = 13.5)
    ret.pheno[[3]] = ret.pheno[[2]]/(ret.pheno[[1]]^2) 
  }  else if(type == "1+2=3" | type == "1+2=3BIG" | type == "1+2=3BIG2")
  {
    #If 3 is a linear combination of 1 and 2
    #Model 1 and 2 as independent
    ret.pheno[[1]] = genotype %*% snp.weights[[1]] 
    ret.pheno[[2]] = genotype %*% snp.weights[[2]]
    ret.pheno[[3]] = 0.5*ret.pheno[[1]] +  0.5*ret.pheno[[2]] 
  }else if(type == "1+2=3,4~5"){
    ret.pheno[[1]] = genotype %*% snp.weights[[1]] 
    ret.pheno[[2]] = genotype %*% snp.weights[[2]]
    ret.pheno[[3]] = 0.5*ret.pheno[[1]] +  0.5*ret.pheno[[2]] 
    ret.pheno[[4]] = genotype %*% snp.weights[[4]]
    ret.pheno[[5]] = genotype %*% snp.weights[[5]]
  }  else if(type == "1,2,3" | type == "1O,2O,3+O") #allows for some overwlapping snps between them.
  {
    #If 3 is a linear combination of 1 and 2
    #Model 1 and 2 as independent
    ret.pheno[[1]] = genotype %*% snp.weights[[1]] 
    ret.pheno[[2]] = genotype %*% snp.weights[[2]] 
    if(type == "1O,2O,3+O")
    {
      ret.pheno[[3]] = genotype %*% snp.weights[[3]] + ret.pheno[[1]]
    } else{
      ret.pheno[[3]] = genotype %*% snp.weights[[3]]
    }
    
  } else if(type == "1,2") #allows for some overwlapping snps between them.
  {
    #If 3 is a linear combination of 1 and 2
    #Model 1 and 2 as independent
    ret.pheno[[1]] = genotype %*% snp.weights[[1]]
    ret.pheno[[2]] = genotype %*% snp.weights[[2]]
  }else
  {
    message("proceeding with null version")
  }
  if(substr(type,1,1) == "0")
  {
    s = str_split(type,pattern = ",")[[1]]
    for(i in 1:length(s))
    {
      if(s[i] == "0")
      {
        ret.pheno[[i]] = rnorm(nrow(genotype),0,env_component[1])
      }else
      {
        message("not yet implemented")
      }
      
    }
  }else #add in the noise and get the R2
  {
    print("DEPRECATED:True R^2 report in the 'r2' argument....")
    r2.mat <- NULL
    r2.mat <- rbind(r2.mat, c(1,2,(cor(ret.pheno[[1]],ret.pheno[[2]])^2)))
    r2.mat <- rbind(r2.mat, c(1,3,(cor(ret.pheno[[1]],ret.pheno[[3]])^2)))
    r2.mat <- rbind(r2.mat, c(2,3,(cor(ret.pheno[[2]],ret.pheno[[3]])^2)))
    ret.pheno$r2 = r2.mat
    
    print("DEPRECTAED:True correlation report in the 'r' argument....")
    r.mat <- NULL
    r.mat <- rbind(r.mat, c(1,2,(cor(ret.pheno[[1]],ret.pheno[[2]]))))
    r.mat <- rbind(r.mat, c(1,3,(cor(ret.pheno[[1]],ret.pheno[[3]]))))
    r.mat <- rbind(r.mat, c(2,3,(cor(ret.pheno[[2]],ret.pheno[[3]]))))
    ret.pheno$r = r.mat
    if(type == "1+2=3" | type == "1+2=3BIG" | type == "1+2=3BIG2" |type == "1+2=3,4~5")
    {
      #CHANGED 4/19
      #Changed this condition- the phenotype is based on phenotypes, not just the genetic piece.
      #This will likely change the simulation results
      message("Using the updated phenotype weights.")
      ret.pheno[[1]] = ret.pheno[[1]] + rnorm(nrow(genotype), 0, sd = env_component[1])
      ret.pheno[[2]] = ret.pheno[[2]] + rnorm(nrow(genotype), 0, sd = env_component[2])
      #ret.pheno[[3]] = 0.5*ret.pheno[[1]] +  0.5*ret.pheno[[2]]  + rnorm(nrow(genotype), 0, sd = env_component[3])
      #modified for quick test april 26
      #ret.pheno[[3]] = 0.2*ret.pheno[[1]] +  0.8*ret.pheno[[2]]  + rnorm(nrow(genotype), 0, sd = env_component[3])
      #modified for quick test april 27, first test went well
      ret.pheno[[3]] = 0.1*ret.pheno[[1]] +  0.9*ret.pheno[[2]]  + rnorm(nrow(genotype), 0, sd = env_component[3])
      message("env_commponents being used:", paste(env_component))
      if(type == "1+2=3,4~5")
      {
        ret.pheno[[4]] = ret.pheno[[4]]+ rnorm(nrow(genotype), 0, sd = env_component[4])
        ret.pheno[[5]] = ret.pheno[[5]] + rnorm(nrow(genotype), 0, sd = env_component[5])
      }
      
    }else if(type == "bmi_like")
    {
      #Don't need to modify the other sones...
      message("BMI-like")
      ret.pheno[[3]]  =  ret.pheno[[3]] + rnorm(nrow(genotype), 0, sd = env_component[1])
      
    }else
    {
      ret.pheno[[3]] = ret.pheno[[3]] + rnorm(nrow(genotype), 0 ,sd = env_component[1])
      ret.pheno[[1]] = ret.pheno[[1]] + rnorm(nrow(genotype), 0, sd = env_component[2])
      ret.pheno[[2]] = ret.pheno[[2]] + rnorm(nrow(genotype), 0, sd = env_component[3])
    }

  }
  return(ret.pheno)
}

genGenotypes <- function(n.snps, mafs)
{
  paternal <- rbinom(n.snps, 1, mafs)
  maternal <- rbinom(n.snps, 1, mafs)
  paternal + maternal
}

#n.snps- how many active snps are we considering?
#n.groups- not used here
#N.tot- how many total sampls are in the cohort?
#genotypes- pass in genotype matrix
#pheno- phenotype combination, i.e. 
#1+2=3- phenotype 3 is a combination of phenotypes 1 and 2
#size: the effect size- a scalar
#active snps- relevant only in the case of 1,2
#overlap snps- shared snps between them.

relatedOverlapPheno <- function(n.snps, n.groups, N.tot = 10000, genotypes, 
                                pheno = "1+2=3", size = 0.5, active_snps = 30, overlap_snps = 25)
{
  print(N.tot)
  genotypes <- genotypes[,1:n.snps]
  #message(interval_size)
  nonoverlapping <- active_snps-overlap_snps
  active.snps <- list()
  if(pheno =="1+2=3" )
  {
    
    active.snps[[1]] <- c(rep(size, 20), rep(0, n.snps-20))
    active.snps[[2]] <- c(rep(0, 20), rep(size, 20), rep(0, n.snps-40))
    active.snps[[3]] <- c(rep(0, n.snps))
  } else if (pheno =="1O,2O,3+O" )
  {
    #First 2 traits have overlap
    #3rd one has independent SNPs
    active.snps[[1]] <- c(rep(size, 50), rep(0, n.snps-50))
    active.snps[[2]] <- c(rep(0, 25), rep(size, 50), rep(0, n.snps-75))
    active.snps[[3]] <- c(rep(0, 101), rep(size, 25), rep(0, n.snps-126))
    
  }else if (pheno =="1+2=3BIG" )
  {
    
    active.snps[[1]] <- c(rep(size, 50), rep(0, n.snps-50))
    active.snps[[2]] <- c(rep(0, 50), rep(size, 100), rep(0, n.snps-150))
    active.snps[[3]] <- c(rep(0, n.snps))
    
  } else if (pheno =="1+2=3BIG2" )
  {
    active.snps[[1]] <- c(rep(0, 300), rep(size, 25), rep(0, n.snps-325))
    active.snps[[2]] <- c(rep(0, 350), rep(size, 25), rep(0, n.snps-375))
    active.snps[[3]] <- c(rep(0, n.snps))
    
  } else if(pheno == "1,2,3") #case with overlap- 2 similar traits
  {
    active.snps[[1]] <- c(rep(size, 40), rep(0, n.snps-40))
    active.snps[[2]] <- c(rep(0, 20), rep(size, 40), rep(0, n.snps-60))
    active.snps[[3]] <- c(rep(0, n.snps))
  } else if(pheno == "1,2") #case with overlap- 2 similar traits
  {
    active.snps[[1]] <- c(rep(size, active_snps), rep(0, n.snps-active_snps))
    active.snps[[2]] <- c(rep(0, nonoverlapping), rep(size, active_snps), rep(0, n.snps-(active_snps + nonoverlapping)))
  } else if(substr(pheno,1,1) == "0") #first element is null
  {
    s = str_split(type,pattern = ",")[[1]]
    for(i in 1:length(s))
    {
      if(s[i] == "0")
      {
        active.snps[[i]] <- c(rep(0, n.snps))
      }else
      {
        message("not yet implemented")
      }
      
    }
  }
  else{
    message("forgot to specify type...")
  }
  phenotypeBuilder(active.snps, genotypes, 0.5,  type = pheno) #this returns a list of all phenotypes for all individuals.
  
}

#Generate GWAS summary statistics across n.snps with cohort overlap of n.overlap.
#Of the total number of samples (N.tot), you can specify a fixed number of samples for each run to ensure consistent power.
#relatedOverlapBETTER(n.snps, i, 3, y1, N.tot = N_tot, 
#                     genotypes=pop1.geno[,1:n.snps], fixed.count = n.fixed)
relatedOverlapBETTER <- function(n.snps, n.overlap, n.groups, y, N.tot = 10000, 
                           genotypes=genotypes.tab.2, fixed.count = 5000, parallel = 0, ret.all.stats = FALSE)
{
  library(parallel)
  if(n.overlap > 0)
  {
    No <- 1:n.overlap
  }else
  {
    No <- NULL
  }
  interval_size = floor((N.tot-n.overlap)/n.groups)
  g <- list()
  t = n.overlap
  if(t >= n.snps)
  {
    t = n.overlap
  }
    #the number of phenotypes must correspond with the number of traits.
    if(fixed.count  == 0)
    {
      for(i in 1:n.groups)
      {
        message(paste0("group " , i))
        if(interval_size > 0)
        {
          cohort_length <- c(No, (t+1):(t+interval_size))
          print(length(cohort_length))
        } else
        {
          cohort_length  <- No
        }
        
        if(ret.all.stats)
        {
          g[[i]] <- lapply(1:n.snps, function(j) summary(lm(scale(y[[i]][cohort_length]) ~ genotypes[cohort_length, j]))$coef)
        }else
        {
          g[[i]] <- sapply(1:n.snps, function(j) summary(lm(scale(y[[i]][cohort_length]) ~ genotypes[cohort_length, j]))$coef[6])
        }
        
        t= t+interval_size
      }
    }else #fixed count implies that each cohort has the exact same number of indivudals
    {
      No <- 1:n.overlap
      if(n.overlap == 0)
      {
        No <- 0
      }
      add <- fixed.count - n.overlap -1 #Need this so we don't overshoot
      
      for(i in 1:n.groups)
      {
        #print("Starting..")
        #print(i)
        start_point = (i-1) * add + n.overlap + (i-1) + as.numeric(n.overlap == 0)
        #print(paste("Start point", start_point))
        cohort_length <- c(No, start_point:(start_point + add))
        #print(paste("last element in list...", cohort_length[fixed.count]))
        #print(cohort_length)
        if(cohort_length[fixed.count] > N.tot)
        {
          message("Error: insufficient number of samples to have this level of overlap.")
        }
        if(n.overlap == 0)
        {
          cohort_length <- cohort_length[-1]
        }
        #print(length(cohort_length))
        if(parallel >= 1)
        {
          parallel.function <- function(j) {
            summary(lm(scale(y[[i]][cohort_length]) ~ genotypes[cohort_length, j]))$coef[6]
        }
          
          g[[i]] <- mclapply( 1:n.snps, FUN=parallel.function, mc.cores = parallel )
          
        }else
        {
          message("Start:", cohort_length[1], "....End: ", cohort_length[length(cohort_length)], "....Length: ", length(cohort_length))
          if(ret.all.stats)
          {
            g[[i]] <- lapply(1:n.snps, function(j) summary(lm(scale(y[[i]][cohort_length]) ~ genotypes[cohort_length, j]))$coef)
          }else
          {
            g[[i]] <- sapply(1:n.snps, function(j) summary(lm(scale(y[[i]][cohort_length]) ~ genotypes[cohort_length, j]))$coef[6])
          }
          
        }
        
      }
    }
  if(ret.all.stats)
  {
    se <- do.call("cbind", lapply(1:length(g), function(i) sapply(g[[i]], function(x) {x[4]})))
    beta <- do.call("cbind", lapply(1:length(g), function(i) sapply(g[[i]], function(x) x[2])))
    #return(do.call("cbind", g))
    return(list("se"=se, "beta"=beta))
  }else
  {
    return(g)
  }
}

####Get covariance data
##input: the phenotype vectors, the overlap counts, and the number of groups
#    first <- relatedOverlapBETTER(n.snps, i, 3, y1, N.tot = N_tot, 
#
#n.snps- how many SNPs overall
#n.overlap -vector of overlap per cohort
#n.groups- how many cohorts we are assessing for 
#n.phenos- how many phenos we are asssessing for, must be the same for all groupsv
#N.tot: a list of how many samples per study by cohort
#y.list- a list of lists, 1 for each group
#y.list <- list(y1,y2,y3)
#n.groups=3
#n.pheno=3
#n.overlap=5000
#N.tot=c()
#5/8- ADDED SCALING TO ALL
calcCovar <- function(n.snps, n.overlap, n.groups,n.pheno=3, y.list, N.tot = 10000, fixed.count = 5000)
{
  
  if(n.overlap > 0)
  {
    No <- 1:n.overlap
  }else
  {
    No <- NULL
    pop.cor  <- lapply(1:n.groups, function(x)  diag(n.pheno))
    pop.se  <- lapply(1:n.groups, function(x)  matrix(0,nrow=n.pheno,ncol=n.pheno))
    return(list("sample.cor"=pop.cor, "sample.se"=pop.se))
  }

  #Loop on all cohorts
  pheno.cor <- list()
  pheno.se <- list()
  pop.se <- list()
  pop.cor <- list()
  #Across all cohorts, for all phenotypes
  for(pop_i in 1:length(y.list))
  {
      #This loop calculates the correlation matrix for the phenotypes (rho_o) and the SE for that estimate
      cohort=y.list[[pop_i]]
      overlapping.phenos <- lapply(cohort[1:n.pheno], function(x) x[No]) #Get just the subset of overlapping samples
      inner.cor.est <- lapply(overlapping.phenos, function(x) lapply(overlapping.phenos, function(y) cor.test(x,y)$estimate)) #these numbers seem off.
      inner.cor.ci <- lapply(overlapping.phenos, function(x) lapply(overlapping.phenos, function(y) {vals <- unlist(cor.test(x,y)$conf.int); (vals[2]-vals[1])/(qnorm(0.975)*2) })) 
      pheno.cor[[pop_i]] <- matrix(as.numeric(do.call("rbind", inner.cor.est)), nrow = n.pheno)
      pheno.se[[pop_i]] <- matrix(as.numeric(do.call("rbind", inner.cor.ci)), nrow = n.pheno)
  }
      #get the final_correlation
      #\frac{No /rho}{sqrt(N_1 N_2)} -->  
    pop.se <- lapply(pheno.se, function(x) (length(No)/sqrt(N.tot*N.tot)) * x)
    pop.cor  <- lapply(pheno.cor, function(x) (length(No)/sqrt(N.tot*N.tot)) * x) %>% lapply(., function(x) {diag(x) <- 1; x})
 
return(list("sample.cor"=pop.cor, "sample.se"=pop.se))
}

fullifyMatrix <- function(m)
{
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

#Statistical helper functions:
getPvals <- function(overlap.list, Z, pcas, npcs =  3)
{
  p.values.overlap <- list()
  for(j in seq_along(overlap.list))
  {
    curr <- NULL
    for(row in 1:nrow(Z[[j]]))
    {
      curr <- rbind(curr, summary(lm(as.matrix(Z[[j]][row,]) ~ as.matrix(pcas[[j]]$v[,1:npcs])))$coef[-1,4])
    }
    p.values.overlap[[j]] <- curr
    print(j)
  }
  return(p.values.overlap)
}

#Statistical helper functions:
getPvalsRow <- function(Z, latent.component)
{
    curr <- NULL  
  for(row in 1:nrow(Z))
    {
      curr <- rbind(curr, summary(lm(as.matrix(Z[row,]) ~ as.matrix(latent.component)))$coef[-1,4])
  }

  return(curr)
}

getPvalsCol <- function(Z, latent.component)
{
  curr <- NULL  
  for(col in 1:ncol(Z))
  {
    curr <- rbind(curr, summary(lm(as.matrix(Z[,col]) ~ as.matrix(latent.component)))$coef[-1,4])
  }
  
  return(curr)
}


getFDRSig <- function(overlap.list, pvals, thresh=0.1, n.traits = 3)
{
  message("Getting BH FDR")
  fdr.list <- list()
  for(j in 1:n.traits){
    print(j)
    fdr.list[[j]] <- list()
    for(i in seq_along(overlap.list))
    {
      fdr <- p.adjust(pvals[[i]][,j], method = "fdr") #in teh gene utils
      fdr.list[[j]][[i]] <- (which(fdr < thresh))
      #print(length(which(p.values.overlap[[i]][,3] < 0.01)))
    } 
  }
  fdr.list
}


##QQPLOT
ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}


##performance metrics
fpr <- function(l, tps, n.snps = 5000) {ns = which(!((1:n.snps) %in% tps)); fps=sum(!(l %in% tps)); fps / (fps + sum(!(ns %in% l)))}
tpr <- function(l, tps) {sum((l %in% tps)) / length(tps)} #TP/TP+FN.
fpp <- function(l) {if(length(l) == 0) {return(0)} else{sum(!(l %in% tps))/length(l)} }

falsePositiveMetrics <- function(fdr.list, tps, n.snps = 5000, N = 15000)
{
  df.fps <- data.frame("fps_2" = sapply(fdr.list[[2]], function(x) sum(!(x %in% tps))),"fps_1" = sapply(fdr.list[[1]], function(x) sum(!(x %in% tps))),"fps_3" = sapply(fdr.list[[3]], function(x) sum(!(x %in% tps)))) %>% mutate("percent_overlap" = overlap.list/N) %>% pivot_longer(cols = c(fps_2, fps_1, fps_3))
  
  df.fpr <- data.frame("fps_2" = sapply(fdr.list[[2]], function(x) fpr(x, tps, n.snps = n.snps)),
                       "fps_1" = sapply(fdr.list[[1]], function(x) fpr(x, tps, n.snps = n.snps)),
                       "fps_3" = sapply(fdr.list[[3]], function(x) fpr(x, tps, n.snps = n.snps))) %>% mutate("percent_overlap" = overlap.list/N) %>% pivot_longer(cols = c(fps_2, fps_1, fps_3))
  
  df.fpp <- data.frame("fps_2" = sapply(fdr.list[[2]], fpp),"fps_1" = sapply(fdr.list[[1]], fpp),"fps_3" = sapply(fdr.list[[3]], fpp)) %>% mutate("percent_overlap" = overlap.list/N) %>% pivot_longer(cols = c(fps_2, fps_1, fps_3))
  return(list("fp_count" = df.fps, "fp_rate" = df.fpr, "fp_proportion" = df.fpp))
  
} 
falsePositiveMetricsBIG <- function(fdr.list,overlaps, tps, N = 15000)
{
  df.fps <- data.frame("fps_2" = sapply(fdr.list[[2]], function(x) sum(!(x %in% tps))),
                       "fps_1" = sapply(fdr.list[[1]], function(x) sum(!(x %in% tps))),
                       "fps_3" = sapply(fdr.list[[3]], function(x) sum(!(x %in% tps)))) %>% 
    mutate("overlap_1" = overlaps[,1], "overlap_2" = overlaps[,2]) %>% pivot_longer(cols = c(fps_2, fps_1, fps_3))
  
  df.fpr <- data.frame("fps_2" = sapply(fdr.list[[2]], function(x) fpr(x, tps)),
                       "fps_1" = sapply(fdr.list[[1]], function(x) fpr(x, tps)),
                       "fps_3" = sapply(fdr.list[[3]], function(x) fpr(x, tps))) %>% 
    mutate("overlap_1" = overlaps[,1], "overlap_2" = overlaps[,2]) %>% pivot_longer(cols = c(fps_2, fps_1, fps_3))
  
  df.fpp <- data.frame("fps_2" = sapply(fdr.list[[2]], fpp),
                       "fps_1" = sapply(fdr.list[[1]], fpp),
                       "fps_3" = sapply(fdr.list[[3]], fpp)) %>% 
    mutate("overlap_1" = overlaps[,1], "overlap_2" = overlaps[,2]) %>% pivot_longer(cols = c(fps_2, fps_1, fps_3))
  return(list("fp_count" = df.fps, "fp_rate" = df.fpr, "fp_proportion" = df.fpp))
  
}
#get  table telling the true positive count at each level (tpc) and the true positive rate at each level (tpr)
truePositiveMetricsBIG <- function(fdr.list,overlaps, tps, N = 15000)
{
  df.tpc <- data.frame("tps_2" = sapply(fdr.list[[2]], function(x) sum((x %in% tps))),
                       "tps_1" = sapply(fdr.list[[1]], function(x) sum((x %in% tps))),
                       "tps_3" = sapply(fdr.list[[3]], function(x) sum((x %in% tps)))) %>% 
    mutate("overlap_1" = overlaps[,1], "overlap_2" = overlaps[,2])  %>% 
    pivot_longer(cols = c(tps_2, tps_1, tps_3))
  
  df.tpr <- data.frame("tps_2" = sapply(fdr.list[[2]], function(x) tpr(x, tps)),
                       "tps_1" = sapply(fdr.list[[1]], function(x) tpr(x, tps)),
                       "tps_3" = sapply(fdr.list[[3]], function(x) tpr(x, tps))) %>%
    mutate("overlap_1" = overlaps[,1], "overlap_2" = overlaps[,2])  %>% pivot_longer(cols = c(tps_2, tps_1, tps_3))
  
  
  return(list("tpr" = df.tpr, "tpc" = df.tpc))
}

truePositiveMetrics <- function(fdr.list, tps, N = 15000)
{
  df.tpc <- data.frame("tps_2" = sapply(fdr.list[[2]], function(x) sum((x %in% tps))),"tps_1" = sapply(fdr.list[[1]], function(x) sum((x %in% tps))),"tps_3" = sapply(fdr.list[[3]], function(x) sum((x %in% tps)))) %>% mutate("percent_overlap" = overlap.list/N) %>% pivot_longer(cols = c(tps_2, tps_1, tps_3))
  
  df.tpr <- data.frame("tps_2" = sapply(fdr.list[[2]], function(x) tpr(x, tps)),"tps_1" = sapply(fdr.list[[1]], function(x) tpr(x, tps)),"tps_3" = sapply(fdr.list[[3]], function(x) tpr(x, tps))) %>% mutate("percent_overlap" = overlap.list/N) %>% pivot_longer(cols = c(tps_2, tps_1, tps_3))
  
  
  return(list("tpr" = df.tpr, "tpc" = df.tpc))
}

lookupCor <- function(str, str2, cor.mat = pheno.cor.mat, multi = FALSE)
{
  f = as.numeric(substr(as.character(str), 4, 4))
  s = as.numeric(substr(as.character(str2), 4, 4))
  if(multi)
  {
    pop1 = as.numeric(substr(as.character(str), 2, 2)) 
    pop2 = as.numeric(substr(as.character(str), 2, 2)) 
    if(pop1 == pop2)
    {
      return(cor.mat[[pop1]][f,s])
    } else
    {
      return(0)
    }
  }else{
    cor.mat[f,s]
  }
 
}

overlapRet <- function(overlap, str, str2, tot = 15000)
{
  f = as.numeric(substr(as.character(str), 2, 2))
  s = as.numeric(substr(as.character(str2), 2, 2))
  if(s == f) #same cohort
  {
    if(as.numeric(substr(as.character(str), 4, 4)) == as.numeric(substr(as.character(str2), 4, 4)))
    {
      return(tot)
    }
    return(overlap)
  }
  return(0)
}


overlapRetBig <- function(overlap, str, str2, tot = 15000)
{
  f = as.numeric(substr(as.character(str), 2, 2))
  s = as.numeric(substr(as.character(str2), 2, 2))
  if(s == f) #same cohort
  {
    if(as.numeric(substr(as.character(str), 4, 4)) == as.numeric(substr(as.character(str2), 4, 4)))
    {
      return(tot)
    }
    return(overlap[s]) #the overlap of whatever cohort they are in.
  }
  return(0)
}
#generate genotypes based on a number of SNPs and existing MAFs
genGenotypes <- function(n.snps, mafs)
{
  paternal <- rbinom(n.snps, 1, mafs)
  maternal <- rbinom(n.snps, 1, mafs)
  paternal + maternal
}



#LDSC related functions have been moved to postMungeUtils.R


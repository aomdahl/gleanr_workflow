library(devtools)
library(optparse)
PSETTINGS <- 10
# Slurm script creation function

create_slurm_script <- function(argv, script_path, startup_commands, additional_args) {
  ###Build the gleaner commadn:
  nargs <- length(argv)
  args.add <- argv[PSETTINGS:nargs]
  arg.names <- paste0("--", names(argv[PSETTINGS:nargs]))
  call.list <- c()
  for(i in 1:length(args.add))
  {
    if( args.add[[i]] != "") #don't add empty stuff.
    {
      if(as.character(args.add[[i]]) %in% c("TRUE", "FALSE"))
      {
        if(as.character(args.add[[i]]) != FALSE)
        {
          
          #print("omitting true/false calls")
          call.list <- c(call.list, paste(arg.names[[i]]))
        }

      }else
      {
        if(arg.names[[i]] == "--outdir")
           {
             args.add[[i]] <- paste0(args.add[[i]], argv$job_name)
        }
      	if(arg.names[[i]] == "--nfactors")
	{
		print(args.add[[i]])
	}
        call.list <- c(call.list, paste(arg.names[[i]], args.add[[i]]))
      }
     
    }
    
  }
  
  slurm_script <- paste0(
    "#!/bin/bash\n",
    "#SBATCH --job-name=", argv$job_name, "\n",
    "#SBATCH --time=",  argv$time, "\n",
    "#SBATCH --partition=",  argv$partition, "\n",
    "#SBATCH --mem=",  argv$memory, "\n",
    "#SBATCH --output=", argv$outdir, "/%j.%x.out\n",
    "#SBATCH --error=",  argv$outdir, "/%j.%x.err\n")
  if(argv$user_account != "")
  {
    #SBATCH -A PI-userid_bigmem+
    slurm_script <- paste0(slurm_script,
                           "#SBATCH -A ", argv$user_account, "\n")
  }
  slurm_script <- paste0(slurm_script,
    "set -e\n",
    "mkdir -p ",argv$outdir,"\n",
    startup_commands,
    "Rscript ", script_path, " ",
    paste(call.list, collapse = " "))

  if(!grepl(slurm_script, pattern="--model_select"))
  {
    slurm_script <- paste0(slurm_script, " --model_select\n")
  }else
  {
    slurm_script <- paste0(slurm_script,"\n")
  }
  
  return(slurm_script)
}

getM <- function(args_input)
{
  setwd("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwasMF/"); devtools::load_all()
  args <-fillDefaultSettings(args_input)
  
  writeRunReport(args)
  
  #Establish the internal settings
  option <- readInSettings(args)
  option$alpha1 <- NA
  option$lambda1 <- NA
  outdir <-args$outdir
  
  #Read in the hyperparameters to explore
  hp <- readInParamterSpace(args)
  #args$WLgamma <- 0
  #args$block_covar <- 0.2 This is the issue
  input.dat <- readInData(args)
  X <- input.dat$X
  return(ncol(X))
}
writeSlurmScript <- function(lines,outdir, job_name )
{
  slurm_file <- paste0(outdir, "/", job_name, ".slurm")
  writeLines(lines, con = slurm_file)
  return(slurm_file)
}
#buildRunScripts(unlist(next.params.to.try$K), args_input)
buildRunScripts <- function(k.iter, args_input){
  init.job.name <- args_input$job_name
  for(K in k.iter)
  {
	args_input$nfactors = K
    print(args_input$nfactors)
    args_input$job_name <- paste0(init.job.name, "_K", K)
    slurm_command <- create_slurm_script(args_input,
                                         script_path = "src/gleaner_run.R",
                                         args_input$init_commands,
                                         additional_args = commandArgs(trailingOnly = TRUE)
                                         
    )
    ofile <- writeSlurmScript(slurm_command,args_input$outdir, args_input$job_name)
    cat("Slurm script created at:", ofile, "\n")
    if(args_input$submit_jobs)
    {
      system(paste0("sbatch ",ofile ))
    }
  }
  args_input$job_name <- init.job.name
}

startup_commands= "source /data/apps/go.sh \n cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/\nml anaconda\nconda activate renv\n"

# Option parsing for command line arguments
option_list <- list(
  make_option(c("--job_name"), type = "character", help = "Job name", default="my_analysis_job"),
  make_option(c("--time"), type = "character", help = "Specify the time for the job in Slurm format, e.g., 01:00:00"),
  make_option(c("--partition"), type = "character", help = "Specify the Slurm partition to use", default = "shared"),
  make_option(c("--memory"), type = "character", help = "Specify the memory for the job, e.g., 4G", default = "4G"),
  make_option(c("--grid_K"), type = "character", help = "Specify which GRID settings to explore, default is to set internally", default = ""),
  make_option(c("--init_commands"), type = "character", help="String with stuff to add to front of bash files",
              default=startup_commands),
  make_option(c("--submit_jobs"), type = "logical", help="Submit the commands after written",
              default=FALSE, action = "store_true"),
  make_option(c("--task"), type = "character", help="tasks to perform, options are BUILD_JOBS, SCORE_JOBS",
              default="BUILD_JOBS"),
  make_option(c("--user_account"), type = "character", default = "", help = "account name for bigmem"),
  make_option(c("--gwas_effects"), type = 'character', help = "Specify the Z or B file, depending on specified weighting scheme. First column is ids of each variant, column names specify the trait"),
  make_option(c("--uncertainty"), type = 'character', help = "Specify the path to the SE or other uncertainty file, depending on the weightin scheme.irst column is ids of each variant, column names specify the trait"),
  make_option(c("--lambda_gc"), type = 'character', help = "Specify the path to the genomic correction coefficients. If none provided, none used", default = ""),
  make_option(c("--trait_names"), type = 'character', help = "Human readable trait names, used for plotting. Ensure order is the same as the order in the input tables.", default = ""),
  make_option(c("--weighting_scheme"), type = 'character', help = "Specify either Z, B, B_SE, B_MAF", default = "B_SE"),
  make_option(c("--covar_matrix"), type = 'character', help = "Path to LDSC estimates of covariance effect. No adjustment made if none provided", default = ""),
  make_option(c("-c", "--converged_obj_change"), type = 'numeric', help = "Specify the objective percent change required to achieve convergence", default = 0.001),
  make_option(c("-i", "--niter"), type = 'numeric', help = "Cap the number of iterations on the matrix learning step", default = 300),
  make_option(c("--outdir"), type = "character", help = "Source file location"),
  make_option(c("--drop_phenos"), type = "character", help = "Specify phenotypes to exclude (useful if encountering issues with covariance adjustment)", default = ""),
  make_option(c("--fixed_first"), type = "logical", help = "if want to remove L1 prior on first factor", action = "store_true", default = FALSE),
  make_option(c("--debug"), type = "logical", help = "if want debug run", action = "store_true", default = FALSE),
  make_option(c("--overview_plots"), type = "logical", help = "To include plots showing the objective, sparsity, etc for each run", action = "store_true", default = FALSE),
  make_option(c("-K", "--nfactors"), type = "character", help = "specify the number of factors. Options are a number, or MAX, KAISER, K-2, CG", default = "MAX"),
  make_option(c("--genomic_correction"), type="character", default= "", help="Specify path to genomic correction data, one per snp.TODO: Also has the flexibility to expand"),
  make_option(c("--bic_var"), type = 'character', help = "Specify the bic method to use. Options are [sklearn_eBIC,sklearn,Zou_eBIC,Zou,dev_eBIC,dev, std,NONE]", default = "sklearn_eBIC"),
  make_option(c("-p", "--param_conv_criteria"), type = 'character', help = "Specify the convergene criteria for parameter selection", default = "BIC.change"),
  make_option(c("-r", "--rg_ref"), type = 'character', help = "Specify a matrix of estimated genetic correlations to initialize V from", default = ""),
  make_option(c("-v", "--verbosity"), type="integer", default= 0, help="How much output information to give in report? 0 is quiet, 1 is loud"),
  make_option(c("-n", "--ncores"), type="integer", default= 1, help="Do you want to run on multiple cores? Only influences the K initialization step at this point."),
  make_option(c("-s", "--sample_sd"), type="character", default= "", help="File containing the standard deviation of SNPs; if provided, used to scale LDSC gcov terms."),
  make_option(c("-b", "--block_covar"), type = "numeric", default= 0.2, help="Specify the degree to which a block structure is enforced, by cluster distance. Default is 0.2"),
  make_option(c("--covar_se_matrix"), type = "character", default="", help="Path to covar se matrix, needed if using Strimmer gamma specification."),
  make_option(c("--intermediate_bic"), type = "character", default="", help="OPTIONAL:Path to intermediate BIC file, will initiate factorization from here. "),
  make_option(c("--subset_seed"), type = "numeric", default=-1, help="Specify a seed to subset the data for a variant test. For debugging purposes."),
  make_option(c("--model_select"), type = "logical", default=FALSE,action ="store_true", help="Specify this if you only wish to do the model selection step."),
  make_option(c("-g", "--WLgamma"), type="character", default= "0", 
              help="Specify the extent of the WL shrinkage on the input covariance matrix.\nCan pass as a number (1 for no covariance adjustment, 0 for full adjustment), or to specify a method (either MLE or Strimmer)")
)


#for debugging purposes
t= c("--gwas_effects=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/finngen_benchmark_2_conservative.beta.tsv",
     "--uncertainty=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//finngen_benchmark_2_conservative.se.tsv",
     "--trait_names=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2//pheno.names.txt",
     "--outdir=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/dev_max_covar_manual", 
     "--fixed_first",  "--nfactors=GRID", "--bic_var=sklearn","--verbosity=1","--ncores=4", "--subset_seed=1000",
     "--covar_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark_2/summary_data/gcov_int.tab.csv",
     "--sample_sd=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/finngen_benchmark_2/sample_sd_report.tsv",
     "--WLgamma=Strimmer", "--covar_se_matrix=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/ldsr_results/finngen_benchmark_2/summary_data/gcov_int_se.tab.csv",
     "--job_name=test", "--time=20", "--partition=shared", "--memory=20G")
fp="/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/"
p=paste0(fp,"gwas_extracts/panUKBB_complete/")
t= c(paste0("--gwas_effects=",p,"/panUKBB_complete_clumped_r2-0.2.beta.tsv"),
     paste0("--uncertainty=",p,"panUKBB_complete_clumped_r2-0.2.se.tsv"),
     paste0("--trait_names=",p,"/panUKBB_complete.trait_list.tsv"),
     paste0("--outdir=",fp, "/results/panUKBB_complete_41K_GRID_dev/"),
     "--fixed_first",  "--nfactors=GRID", "--bic_var=dev",
     paste0("--covar_matrix=", fp, "/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv"),
     paste0("--sample_sd=",p,"sample_sd_report.tsv"),"--WLgamma=Strimmer", 
      paste0("--covar_se_matrix=",fp,"/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv"),
     "--job_name=dev-grid-search", "--time=24:00:00", "--partition=shared", "--memory=75G")

p=paste0(fp,"gwas_extracts/panUKBB_complete/")
t= c(paste0("--gwas_effects=",p,"/panUKBB_complete_clumped_r2-0.2.beta.tsv"),
     paste0("--uncertainty=",p,"panUKBB_complete_clumped_r2-0.2.se.tsv"),
     paste0("--trait_names=",p,"/panUKBB_complete.trait_list.tsv"),
     paste0("--outdir=", fp, "/results/panUKBB_complete_41K_sklearn_eBIC/"),
     "--fixed_first",  "--nfactors=GRID", "--bic_var=sklearn_eBI",
     paste0("--covar_matrix=", fp, "/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv"),
     paste0("--sample_sd=",p,"sample_sd_report.tsv"),"--WLgamma=Strimmer", 
     paste0("--covar_se_matrix=",fp,"/ldsr_results/panUKBB_complete/summary_data/gcov_int_se.tab.csv"),
     "--job_name=dev-grid-search", "--time=24:00:00", "--partition=shared", "--memory=75G")




#args_input <- parse_args(OptionParser(option_list = option_list),args = t)
args_input <- parse_args(OptionParser(option_list = option_list))

additional_args <- commandArgs(trailingOnly = TRUE)

if(args_input$grid_K == "")
{
  M <- suppressMessages(getM(args_input))
  k.range <- 1:(M-1)
}


if(args_input$task == "BUILD_JOBS")
{
  init.range=4
  if(args_input$grid_K == "")
  {
    k.iter <- unlist(list("K" = floor(quantile(k.range,seq(0, 1, 1/init.range)))[-1]))
  }else
  {
    k.iter <- unlist(as.numeric(unlist(strsplit(args_input$grid_K,"," ))))
  }
  # Call the function to create the Slurm script
 
  buildRunScripts(k.iter, args_input)
  
}
if(grepl(args_input$task, "SCORE"))
{
  files <- list.files(path = args_input$outdir, pattern = "*_bic_dat.RData")
  valid.k.measures <- which(grepl(files, pattern = "K[0-9]+"))
  performance.dat <- NULL
  ret.list <- list()
  for(f in files[valid.k.measures])
  {
    k.val <- stringr::str_match(f, pattern = "K([0-9]+)_")[2]
    load(paste0(args_input$outdir,f))
    ret.list[[f]] <- bic.dat
    performance.dat <- rbind(performance.dat, 
                             c(k.val, bic.dat$min.dat$min_sum, bic.dat$min.dat$alpha, bic.dat$min.dat$lambda))
  }
  performance.dat <- data.frame(performance.dat) %>% magrittr::set_colnames(c("Kinit", "BIC_sum","alpha", "lambda")) %>% mutate_all(as.numeric)
  
  performance.dat$bic_global <- getGlobalBICSum(ret.list)
  grid_record <- list("curr_runs"=performance.dat, "test_dat"=ret.list)
  next.params.dat <- gatherSearchData(ret.list,k.range,grid_record)
  next.params.to.try <- next.params.dat$next_params
  #best_k = performance.dat$K[which.min(performance.dat$bics_global)]
  #grid <- list("BIC"=performance.dat$BIC_sum, "K"=performance.dat$K)
  #next.params.to.try <- chooseNextParams(best_k, grid,k.range)
  
  if(grepl(args_input$task, "BUILD"))
  {
    buildRunScripts(unlist(next.params.to.try$K), args_input)
  }

}

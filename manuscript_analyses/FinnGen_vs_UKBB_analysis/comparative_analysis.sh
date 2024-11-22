#!/bin/bash

#######################################################################
#
#  Script and commands for performing cross-cohort analysis on GWAS from the UKBB and FinnGen reported in manuscript
#
#	Ashton Omdahl, January 2024
# 	These commandes were run using on an HPC `slurm` workload manager and so utilize functions like `sbatch`
#	Analysis was performed on the Johns Hopkins ARCH Rockfish computing cluster (https://www.arch.jhu.edu/about-arch/)
#   Script step summary:
#		- Steps 1-7: Downloading and extracting GWAS summary statistics data
#		- Steps 8-11: Performing XT-LDSC on 194 traits to manually filter out those with pairwise rg > 0.7
#		- Steps 12-15: Clumping SNPs from selected traits and extracting summary statistics as matrices
#		- Steps 16-17: GLEANR analysis
#		- Steps 18-: 
#######################################################################

#File download info from a local file: ukbb_fingne_phenotypes.Rmd
# Commands to download files generated there

#0) Reference information:
RUN_SCRIPTS=run_scripts/finngen_ukbb_benchmark/v2_expanded_run/covar_influence/
FG_DATA=/gwas_extracts/finngen_benchmark_2
UK_DATA=/gwas_extracts/ukbb_benchmark_2
ONEKG_REF=/scratch16/abattle4/ashton/prs_dev/1000genomes_refLD/
RUN_DIR=manuscript_analyses/FinnGen_vs_UKBB_analysis
Starts with downlaoded....
Stats with ukbb_fingen_phenotypes.Rmd for selection...

#After manually + computationally revewing manifests, specified studies in both FinnGen and UKBB for consideration

#1) Download the files:
   bash ${RUN_DIR}/finngen_download_commands_jan2024.sh
   bash ${RUN_DIR}/ukbb_saige_download_commands_jan2024.sh

#2) Manually construct lists for analysis
        #Format source GWAS into an input file with a tab-delimited file with columns of PATH TRAIT_NAME GENOME_BUILD #CASE/#CONTROL COHORT
        #These were stored in trait_selections/finngen_benchmark_2.studies.tsv and trait_selections/ukbb_benchmark_2.studies.tsv,
        # and are later moved into trait_selections/finngen_benchmark_2.FULL_OLD.studies.tsv and trait_selections/ukbb_benchmark_2.FULL_OLD.studies.tsv, respectively

#3) Format GWAS for XT-LDSR analysis
        snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/finngen_comparative_2_config.yaml -j 4 gwas_extracts/finngen_benchmark_2/missingness_report.tsv  --profile profiles/rockfish/

        snakemake --snakefile rules/ldsr_pairwise.smk --configfile config/ukbb_comparative_2_config.yaml -j 4 gwas_extracts/ukbb_benchmark_2/missingness_report.tsv --dry-run  --profile profiles/rockfish/

#4) Run XT-LDSR on matched phenotypes from each cohort to select only those that are strongly correlated (rg > 0.8)
        #A) Build the reference files
                UK=gwas_extracts/ukbb_benchmark_2/munged_gwas_file_paths.txt
                FG=gwas_extracts/finngen_benchmark_2/munged_gwas_file_paths.txt
                realpath gwas_extracts/ukbb_benchmark_2/*.sumstats.gz > $UK
                realpath gwas_extracts/finngen_benchmark_2/*.sumstats.gz > $FG
                OUT=/scratch16/abattle4/ashton/snp_networks/scratch/finngen_v_ukbb/v2_expanded_run/rg_cross_cohort/
        #B) Build scripts to run on :
                bash src/xt-ldsr_across_two_lists.sh $OUT $UK $FG
        #C) Run each cross-trait output.
                for i in {1..79}; do
                        echo $i
                        bash pheno_$i.sh
                done
        #D) Carefully evaluate the output and select which phenotypes should stay, and which should go.
                cd /scratch16/abattle4/ashton/snp_networks/scratch/finngen_v_ukbb/v2_expanded_run/rg_cross_cohort
                grep -A 3 "Summary of Genetic" type_1_diabetes_w_71.log | sed '/^$/d'| tail -n +2| head -n 1 > joined.rg.tsv
                for i in *.log; do
                echo $i
                grep -A 3 "Summary of Genetic" $i | tail -n 2 | sed '/^$/d' >> joined.rg.tsv
                done
                #This will be considered by evaluating: genetic correlation between the phenotypes, significance of those estimates, Neff.
        #E) Evaluate the results, done on my local machine (ukbb_fingen_phenotypes.Rmd) and set thresholds (rg > 0.8, p < 1e-4)
#### followed by some Rmd script on local machine

##########
```{bash}

awk '(FNR == NR) {arr[$1]=$1;next} !($2 in arr) {print $0}' drop.names.txt ukbb_benchmark_2.FULL_OLD.studies.tsv > ukbb_benchmark_2.studies.tsv

awk '(FNR == NR) {arr[$1]=$1;next} !($2 in arr) {print $0}' drop.names.txt finngen_benchmark_2.FULL_OLD.studies.tsv  > finngen_benchmark_2.studies.tsv

```



need to make sure all things in the right order:

```{bash}
UK=ukbb_benchmark_2
FG=finngen_benchmark_2
# mungeCommands
touch gwas_extracts/${UK}/munge_sumstats.all.sh
touch gwas_extracts/${FG}/munge_sumstats.all.sh

#splitMungeCommands
touch gwas_extracts/${UK}/munge_calls/ldsc_line_runs.*.sh
touch gwas_extracts/${FG}/munge_calls/ldsc_line_runs.*.sh

#ldscMunge

touch gwas_extracts/${UK}/*.sumstats.gz
touch gwas_extracts/${FG}/*.sumstats.gz
```

Then run the snakemaker:
```{bash}

snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/finngen_comparative_2_config.yaml -j 1 gwas_extracts/finngen_benchmark_2/missingness_report.tsv

snakemake --snakefile rules/ldsr_pairwise.smk --configfile config/ukbb_comparative_2_config.yaml -j 4 gwas_extracts/ukbb_benchmark_2/missingness_report.tsv


#that worked.

snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/finngen_comparative_2_config.yaml -j 5 --profile profiles/rockfish/  ldsr_results/finngen_benchmark_2/summary_data/gcov_int.tab.csv



snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/ukbb_comparative_2_config.yaml -j 5 --profile profiles/rockfish/  ldsr_results/ukbb_benchmark_2/summary_data/gcov_int.tab.csv

#single runs for debugging purposes
snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/finngen_comparative_2_config.yaml -j 1 ldsr_results/finngen_benchmark_2/summary_data/gcov_int.tab.csv
```
##########

# Run XT-LDSC
snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/finngen_comparative_2_config.yaml -j 5 --profile profiles/rockfish/  ldsr_results/finngen_benchmark_2/summary_data/gcov_int.tab.csv
snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/ukbb_comparative_2_config.yaml -j 5 --profile profiles/rockfish/  ldsr_results/ukbb_benchmark_2/summary_data/gcov_int.tab.csv

#2) Extract the SNPs from these studies for factorization:
#Perform XT-LDSC pariwsie across input traits to estimate covariances structure
snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/finngen_comparative_2_config.yaml -j 5 --profile profiles/rockfish/  ldsr_results/finngen_benchmark_2/summary_data/gcov_int.tab.csv
snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/ukbb_comparative_2_config.yaml -j 5 --profile profiles/rockfish/  ldsr_results/ukbb_benchmark_2/summary_data/gcov_int.tab.csv

# Select out variants for consideration in analysis
snakemake --snakefile rules/extract_factorize.smk  --configfile config/ukbb_comparative_2_config.yaml -j 2   gwas_extracts/ukbb_benchmark_2/ukbb_benchmark_2.1000G.txt --dry-run

snakemake --snakefile rules/extract_factorize.smk  --configfile config/finngen_comparative_2_config.yaml -j 2   gwas_extracts/finngen_benchmark_2/finngen_benchmark_2.1000G.txt 

#.......

#Select out SNPs nominated by both UKBB and FinnGen both for downstream analysis:
    #Done interactively in the notebook: 
    #/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/manuscript_analyses/FinnGen_vs_UKBB_analysis/comparing_raw_SNP_data.Rmd
    [ ]     #TODO ^ Clean this up to just reflect the core analysis we need
    cp ${UK_DATA}/snps_in_both.txt  ${UK_DATA}/ukbb_benchmark_2_conservative.1000G.txt
    cp ${FG_DATA}/snps_in_both.txt  ${FG_DATA}/finngen_benchmark_2_conservative.1000G.txt      

#Prune based on these lists, and then extract the corresponding GWAS information (rules extract_sumstats,/src)

snakemake --snakefile rules/extract_factorize.smk  --configfile config/ukbb_comparative_2_config.yaml -j 5 --profile profiles/rockfish/   gwas_extracts/ukbb_benchmark_2/ukbb_benchmark_2.beta.tsv

snakemake --snakefile rules/extract_factorize.smk  --configfile config/finngen_comparative_2_config.yaml -j 5 --profile profiles/rockfish/   gwas_extracts/finngen_benchmark_2/finngen_benchmark_2.beta.tsv

# Perform GLEANR factorization:
    #Get the clean trait names:
    awk '{print $2}'  trait_selections/ukbb_benchmark_2.studies.tsv > gwas_extracts/finngen_benchmark_2/pheno.names.txt
    cp gwas_extracts/finngen_benchmark_2/pheno.names.txt gwas_extracts/ukbb_benchmark_2/pheno.names.txt

    # With covariance adjustment
    sbatch ${RUN_SCRIPTS}/ukbb_benchmark_2_GRID_sklearn_eBIC.14K.sh
    sbatch ${RUN_SCRIPTS}/finngen_benchmark_2_GRID_sklearn_eBIC.41K.sh
    # Without covariance adjustment
    sbatch ${RUN_SCRIPTS}/finngen_benchmark_2_GRID_sklearn_eBIC.no_adj.41K.sh
    sbatch ${RUN_SCRIPTS}/ukbb_benchmark_2_GRID_sklearn_eBIC.no_adj.41K.sh

#Files 
/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/finngen_benchmark_2/conservative_1e-5/covar_influence
/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/ukbb_benchmark_2/conservative_1e-5/covar_influence/

*BIC-sklearn_eBIC_K-41_no_covarlatent.factors.txt
*BIC-sklearn_eBIC_K-41*
#Compare performance with and without adjustment
fig_3_overlap_and_bar.R

#add in all config files:
config/ukbb_comparative_2_config.yaml
munging scripts
reference data where can yo get it?

ldscref=/data/abattle4/aomdahl1/reference_data/eur_w_ld_chr #Location is hard-coded into scripts....

/data/abattle4/aomdahl1/reference_data/ldsc_reference/1000G_EUR_Phase3_baseline
/data/abattle4/aomdahl1/reference_data/eur_w_ld_chr

/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/run_scripts/finngen_ukbb_benchmark/v2_expanded_run/covar_influence

${RUN_SCRIPTS}/ukbb_benchmark_2_GRID_sklearn_eBIC.14K.sh

files in gwas_decomp_ldsc/src/

/variant_lookup.sh
/snp_cleanup.sh
/variant_to_rsid.sh
update the paths on the thing.

finngen_benchmark_2.FULL_OLD.studies.tsv and 

ukbb_benchmark_2.FULL_OLD.studies.tsv
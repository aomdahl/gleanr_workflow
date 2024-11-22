#  Cross-cohort analysis on GWAS from the UKBB and FinnGen
*Ashton Omdahl, Analysis reported began January 2024, updated Nov. 2024*

These commands were run using on an HPC `slurm` workload manager and so utilize functions like `sbatch`. Analysis was performed on the Johns Hopkins ARCH Rockfish computing cluster (https://www.arch.jhu.edu/about-arch/)
## Script step summary:
- **Steps 1-3**: Downloading and extracting GWAS summary statistics data
- **Steps 4**: Performing XT-LDSR on traits to compare cross-cohort correlation and select phenotypes genetically similar
- **Steps 5-7**: Perform XT-LDSR on paris of traits within a cohort to estimate inflation due to sample sharing (for `C` matrix)
- **Steps 8-9**: Extract SNPs for evaluation in GLEANR
- **Steps 10-11**: GLEANR analysis and evaluation of factor similarity


## Analysis steps
**0) Set reference variables:**
```bash
RUN_SCRIPTS=run_scripts/finngen_ukbb_benchmark/v2_expanded_run/covar_influence/
FG_DATA=/gwas_extracts/finngen_benchmark_2
UK_DATA=/gwas_extracts/ukbb_benchmark_2
ONEKG_REF=/scratch16/abattle4/ashton/prs_dev/1000genomes_refLD/
RUN_DIR=manuscript_analyses/FinnGen_vs_UKBB_analysis
```
**1) Download GWAS files:**
Traits for download were selected based on biobank file manifests through a combination of manually and computational comparisons. This processes is outlined in `manuscript_analyses/FinnGen_vs_UKBB_analysis/analyses/ukbb_finngen_phenotypes_matching.Rmd`
```bash
   bash ${RUN_DIR}/finngen_download_commands_jan2024.sh
   bash ${RUN_DIR}/ukbb_saige_download_commands_jan2024.sh
```
**2) Manually construct lists for analysis**

Here, I formatted source GWAS into an input file with a tab-delimited file with columns of PATH TRAIT_NAME GENOME_BUILD #CASE/#CONTROL COHORT
These were stored in `trait_selections/finngen_benchmark_2.studies.tsv` and `trait_selections/ukbb_benchmark_2.studies.tsv`, and are later moved into `trait_selections/finngen_benchmark_2.FULL_OLD.studies.tsv` and `trait_selections/ukbb_benchmark_2.FULL_OLD.studies.tsv`, respectively

**3) Format GWAS for XT-LDSR analysis**
```bash
ml snakemake
snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/finngen_comparative_2_config.yaml -j 4 gwas_extracts/finngen_benchmark_2/missingness_report.tsv  --profile profiles/rockfish/

snakemake --snakefile rules/ldsr_pairwise.smk --configfile config/ukbb_comparative_2_config.yaml -j 4 gwas_extracts/ukbb_benchmark_2/missingness_report.tsv --dry-run  --profile profiles/rockfish/
```
This command writes LDSR-friendly summary statistics files into `gwas_extracts/finngen_benchmark_2/*.sumstats.gz` and `gwas_extracts/ukbb_benchmark_2/*.sumstats.gz`

**4) Run XT-LDSR on matched phenotypes from each cohort**

The aim here is to select for GLEANR analysis only those that are strongly correlated (e.g. rg > 0.8)

> *4A) Build the reference files, listing all the summary statistics for analysis*

```bash
UK=gwas_extracts/ukbb_benchmark_2/munged_gwas_file_paths.txt
FG=gwas_extracts/finngen_benchmark_2/munged_gwas_file_paths.txt
realpath gwas_extracts/ukbb_benchmark_2/*.sumstats.gz > $UK
realpath gwas_extracts/finngen_benchmark_2/*.sumstats.gz > $FG
OUT=/scratch16/abattle4/ashton/snp_networks/scratch/finngen_v_ukbb/v2_expanded_run/rg_cross_cohort/
```
> 4B) Build the reference files, listing all the summary statistics for analysis
```bash
bash src/xt-ldsr_across_two_lists.sh $OUT $UK $FG
```
> 4C) Run each cross-trait output.
```bash 
for i in {1..79}; do
        echo $i
        bash pheno_$i.sh
done
```
> 4D) Gather output data to evaluate which phenotypes should stay, and which should go.
```bash 
cd $OUT
grep -A 3 "Summary of Genetic" type_1_diabetes_w_71.log | sed '/^$/d'| tail -n +2| head -n 1 > joined.rg.tsv
for i in *.log; do
echo $i
grep -A 3 "Summary of Genetic" $i | tail -n 2 | sed '/^$/d' >> joined.rg.tsv
done
```
> 4E) Evaluate the results, pick final candidates

This was done in the interactive R-notebook, `manuscript_analyses/FinnGen_vs_UKBB_analysis/analyses/ukbb_finngen_xt-ldsr_selection.Rmd`. Traits to omit were placed in the file `drop.names.txt`

**5) Update the trait lists for GLEANR analysis pipeline**
```{bash}
mv trait_selections/ukbb_benchmark_2.studies.tsv trait_selections/ukbb_benchmark_2.FULL_OLD.studies.tsv;
mv trait_selections/finngen_benchmark_2.studies.tsv trait_selections/finngen_benchmark_2.FULL_OLD.studies.tsv;
awk '(FNR == NR) {arr[$1]=$1;next} !($2 in arr) {print $0}' trait_selections/drop.names.txt trait_selections/ukbb_benchmark_2.FULL_OLD.studies.tsv > trait_selections/ukbb_benchmark_2.studies.tsv
awk '(FNR == NR) {arr[$1]=$1;next} !($2 in arr) {print $0}' trait_selections/drop.names.txt trait_selections/finngen_benchmark_2.FULL_OLD.studies.tsv  > trait_selections/finngen_benchmark_2.studies.tsv
#Remove intermediate files so snakemake reruns
rm  gwas_extracts/ukbb_benchmark_2/missingness_report.tsv
rm  gwas_extracts/finngen_benchmark_2/missingness_report.tsv
```

**6) Get Snakemake to rerun QC on just the selected summary statistics**
Requires some updating of intermediate files so we don't re-do extra work
```bash
UK=ukbb_benchmark_2
FG=finngen_benchmark_2
# Rule mungeCommands
touch gwas_extracts/${UK}/munge_sumstats.all.sh; touch gwas_extracts/${FG}/munge_sumstats.all.sh
#Rule splitMungeCommands
touch gwas_extracts/${UK}/munge_calls/ldsc_line_runs.*.sh; touch gwas_extracts/${FG}/munge_calls/ldsc_line_runs.*.sh
#Rule ldscMunge
touch gwas_extracts/${UK}/*.sumstats.gz;touch gwas_extracts/${FG}/*.sumstats.gz
```
Then run the snakemaker:
```bash
snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/finngen_comparative_2_config.yaml -j 1 gwas_extracts/finngen_benchmark_2/missingness_report.tsv

snakemake --snakefile rules/ldsr_pairwise.smk --configfile config/ukbb_comparative_2_config.yaml -j 4 gwas_extracts/ukbb_benchmark_2/missingness_report.tsv
```
**7) Perform XT-LDSR pairwise within cohort to generate cohort overlap estimates (C)**
```bash
snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/finngen_comparative_2_config.yaml -j 5 --profile profiles/rockfish/  ldsr_results/finngen_benchmark_2/summary_data/gcov_int.tab.csv

snakemake --snakefile rules/ldsr_pairwise.smk  --configfile config/ukbb_comparative_2_config.yaml -j 5 --profile profiles/rockfish/  ldsr_results/ukbb_benchmark_2/summary_data/gcov_int.tab.csv
```
**8) Select out candidate variants from GWAS for consideration in analysis (p < 1e-5)**
```bash
snakemake --snakefile rules/extract_factorize.smk  --configfile config/ukbb_comparative_2_config.yaml -j 2   gwas_extracts/ukbb_benchmark_2/ukbb_benchmark_2.1000G.txt

snakemake --snakefile rules/extract_factorize.smk  --configfile config/finngen_comparative_2_config.yaml -j 2   gwas_extracts/finngen_benchmark_2/finngen_benchmark_2.1000G.txt 
```
**9) Filter down to just SNPs nominated by both cohorts for downstream analysis**
This was done interactively in the notebook: `manuscript_analyses/FinnGen_vs_UKBB_analysis/analyses/comparing_raw_SNP_data.Rmd`
    [ ]     #TODO ^ Clean this up to just reflect the core analysis we need
Update the intermediate files to reflect these SNPs:
```bash
    cp ${UK_DATA}/snps_in_both.txt  ${UK_DATA}/ukbb_benchmark_2_conservative.1000G.txt
    cp ${FG_DATA}/snps_in_both.txt  ${FG_DATA}/finngen_benchmark_2_conservative.1000G.txt      
```
**9) Prune based on these lists, and then extract the corresponding GWAS information**
```bash
snakemake --snakefile rules/extract_factorize.smk  --configfile config/ukbb_comparative_2_config.yaml -j 5 --profile profiles/rockfish/   gwas_extracts/ukbb_benchmark_2/ukbb_benchmark_2.beta.tsv
snakemake --snakefile rules/extract_factorize.smk  --configfile config/finngen_comparative_2_config.yaml -j 5 --profile profiles/rockfish/   gwas_extracts/finngen_benchmark_2/finngen_benchmark_2.beta.tsv
```
**10)  Perform GLEANR factorization:**
Get the clean trait names:
```bash
awk '{print $2}'  trait_selections/ukbb_benchmark_2.studies.tsv > gwas_extracts/finngen_benchmark_2/pheno.names.txt
cp gwas_extracts/finngen_benchmark_2/pheno.names.txt gwas_extracts/ukbb_benchmark_2/pheno.names.txt
```    
Run both cohorts *with* covariance adjustment
```
sbatch ${RUN_SCRIPTS}/ukbb_benchmark_2_GRID_sklearn_eBIC.14K.sh
sbatch ${RUN_SCRIPTS}/finngen_benchmark_2_GRID_sklearn_eBIC.41K.sh
```
Run both cohorts *without* covariance adjustment
```bash
sbatch ${RUN_SCRIPTS}/finngen_benchmark_2_GRID_sklearn_eBIC.no_adj.41K.sh
sbatch ${RUN_SCRIPTS}/ukbb_benchmark_2_GRID_sklearn_eBIC.no_adj.41K.sh
```

**11) Evaluate and compare output factorizations:**
Analysis reported in manuscript carried out in `/manuscript_analyses/figure_scripts/fig_3_overlap_and_bar.R`


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

./trait_selections/drop.names.txt
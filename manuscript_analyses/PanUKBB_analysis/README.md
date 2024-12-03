# Analysis of 137 diverse Pan-UKBB GWAS traits
*Ashton Omdahl, January 2024*

These commands were run using on an HPC `slurm` workload manager and so utilize functions like `sbatch`
Analysis was performed on the Johns Hopkins ARCH Rockfish computing cluster (https://www.arch.jhu.edu/about-arch/)
## Script step summary:
- Steps 0-4: Downloading and extracting GWAS summary statistics data
- Steps 5-8: Performing XT-LDSC on 194 traits to manually filter out those with pairwise rg > 0.7
- Steps 9-11: Clumping SNPs from selected traits and extracting summary statistics as matrices
- Steps 12-14: GLEANR analysis
- Steps 15-18: Evaluation and visualization of GLEANR factors

## Analysis steps
**0) Set reference variables:**
```bash
BASE_DIR=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/
PU_SRC=manuscript_analyses/PanUKBB_analysis/src/
GWAS_DIR=/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB
DWNLD_SCRIPTS=${PU_SRC}/download_and_filter_scripts/
EXTRACT_SCRIPTS=${PU_SRC}/extract_summary_stats_scripts/
NOTEBOOKS=manuscript_analyses/PanUKBB_analysis/
SCRATCH=/scratch16/abattle4/ashton/snp_networks/scratch/panUKBB_analysis/
```
**1) Select traits for download, build download commands:**

File download info in `/manuscript_analyses/PanUKBB_analysis/analysis/trait_selection.Rmd`


**2) Download target GWAS files:**

First get the manifests and do your filtering. 

These were interactively downloaded as `.csv` files from `https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288` on Feb 8, 2024.

These were processed in `/manuscript_analyses/PanUKBB_analysis/analysis/trait_selection.Rmd`
```bash
	cd ${GWAS_DIR}
  	sbatch ${DWNLD_SCRIPTS}/download_biomarkers_flat.sh
  	sbatch ${DWNLD_SCRIPTS}/download_phenos_flat.sh
```
**2) Select out the quality variants for consideration downstream**
```bash
	cd ${GWAS_DIR}
	ml anaconda; conda activate renv
	Rscript ${DWNLD_SCRIPTS}/quick_filter.R
```

**3) Extract the European-only summary statistics:**
```bash
  bash ${EXTRACT_SCRIPTS}/build_extract_scripts.sh #This writes out scripts to extract the European sample data per file. Paths are hard-coded in.
  #Launch the commands on slurm		
  bash ${EXTRACT_SCRIPTS}/submit_extractions.sh
```

**4) Manually compile the trait manifest for the 194 traits pipeline.**

This file is a TSV with no header, with columns "PATH TRAIT GENOME_BUILD N COHORT".
File in `trait_selections/panUKBB_complete.studies.tsv`
  
**5) Format GWAS for LDSC analysis**

Launch the premunging script to get these into a format that the LDSC munger will recognize.
```bash
	sbatch manuscript_analyses/PanUKBB_analysis/src/premunge_all.sh
	cd ${BASE_DIR}
	#Clean up headers, extract the right columns
	sed -i 's/NEGLOG10_PVAL_EUR/PVAL/g' gwas_extracts/panUKBB_complete/munge_sumstats.all.sh
	#It didn't add the signed sum stats,nor did it add SE. look into this. But for now, do it manually
	sed -i 's/--keep-maf/--keep-maf --signed-sumstats beta_EUR,0/g'  gwas_extracts/panUKBB_complete/munge_sumstats.all.sh
	sed -i 's/beta_EUR,0/beta_EUR,0 --se se_EUR/g'  gwas_extracts/panUKBB_complete/munge_sumstats.all.sh
	#NOTE: the SE call only works with my moded sumstat munger script.
```
Note that I don't use the generic snakemake pipeline here for munging because these files have an unusual format and I am only running part of the pipeline here, so had to it all manually. Otherwise, could have called it with something like:
```bash
	snakemake --snakefile rules/ldsr_pairwise.smk -j 1 gwas_extracts/panUKBB_complete/sample_sd_report.tsv  --configfile config/panUKBB_complete_config.yaml
```
**6) Reformat GWAS to prepare for XT-LDSC**
First split the munge commands across individual jobs I can submit separately:
```bash
	ml snakemake
	snakemake --snakefile  rules/ldsr_pairwise.smk  -j 1 gwas_extracts/panUKBB_complete/munge_calls/ldsc_line_runs.mean_platelet_volume.sh --configfile config/panUKBB_complete_config.yaml
```

Then submit the munging jobs to run on slurm job manager:
```bash
	snakemake --snakefile rules/ldsr_pairwise.smk gwas_extracts/panUKBB_complete/missingness_report.tsv  -j 5 --configfile config/panUKBB_complete_config.yaml --profile  profiles/rockfish/
```
**7)Run XT-LDSC on GWAS to estimate genetic correlation and correlation due to sample sharing**

Generate the pairwise LDSC commands and run them
```bash
	#cleanup previous runs:
	rm ldsr_results/panUKBB_complete/rg_ldsr/*.sh
	snakemake --snakefile rules/ldsr_pairwise.smk -j 1 ldsr_results/panUKBB_complete/rg_ldsr/hdl_cholesterol_ldsc.run.sh  --configfile config/panUKBB_complete_config.yaml
```
Modify ldsc scripts to work with the format of these files, LDSC output formatting was irregular and we need to update to EUR.:
```bash
	for i in  ` ls ldsr_results/panUKBB_complete/rg_ldsr/*.sh`; do

		sed -i 's/-chr//g' $i
		sed -i 's/UKBB.ALL.ldscore\//UKBB.ALL.ldscore\/UKBB.EUR/g' $i
	done
```
Run the commands:
```bash
snakemake --snakefile rules/ldsr_pairwise.smk -j 5 ldsr_results/panUKBB_complete//summary_data/gcov_int.tab.csv  --configfile config/panUKBB_complete_config.yaml --profile profiles/rockfish
```

**8) Manually evaluate the results, select an rg threshold to filter "redundant" studies**

This was done in an R notebook interactively- `${NOTEBOOKS}/trait_selection_part2_postLDSC.Rmd`

**9) Manually select those files out for downstream consideration**

```bash
cp ${SCRATCH}/rg_filtered_0.7_traits.txt  ldsr_results/panUKBB_complete/rg_filtered_0.7_traits.txt 
awk -F "." '(FNR == NR) {arr[$1];next} ($4 in arr) {print $0}' ldsr_results/panUKBB_complete/rg_filtered_0.7_traits.txt  gwas_extracts/panUKBB_complete/missingness_report.tsv > gwas_extracts/panUKBB_complete/missingness_report.SUB.tsv
mv gwas_extracts/panUKBB_complete/missingness_report.tsv gwas_extracts/panUKBB_complete/missingness_report.COMPLETE.tsv
mv gwas_extracts/panUKBB_complete/missingness_report.SUB.tsv gwas_extracts/panUKBB_complete/missingness_report.tsv
```

**10) Extract variants at target p-value threshold, and calculate pleiotropy scores for clumping procedure**
```bash
ml anaconda
conda activate std
python src/unionVariants.py --gwas_list gwas_extracts/panUKBB_complete/missingness_report.tsv  --output gwas_extracts/panUKBB_complete/panUKBB_complete.union.txt --type ldsc_custom --pval 1e-5 --gwas_dir gwas_extracts/panUKBB_complete/ --output_counts gwas_extracts/panUKBB_complete/panUKBB_complete.union_freq.txt
conda deactivate
cd gwas_extracts/panUKBB_complete/
```
Map the RSID identifiers to chr:pos identifiers for clumping
```bash
GE=/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/gwas_extracts/panUKBB_complete/
cd /scratch16/abattle4/ashton/snp_networks/custom_l1_factorization
ml anaconda
conda activate renv
Rscript ${PU_SRC}/snp_filtering_pre_pruning.R ${GE}/panUKBB_complete.union_freq.txt ${GE}/panUKBB_complete.1000G_freqs.txt
conda deactivate
```

**11) Clump using the UKBB reference and plink**

Run the clumping:
```bash
	ml plink/1.90b6.4
	plink --bfile /scratch16/abattle4/ashton/prs_dev/1000genomes_refLD/ref --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 250 --clump  gwas_extracts/panUKBB_complete/panUKBB_complete.1000G_freqs.txt --out gwas_extracts/panUKBB_complete/panUKBB_complete_clumping_250kb_r2_0.2
```
Convert it into a format that works with the snakemake pipeline
```bash
	awk '{print $3}'   gwas_extracts/panUKBB_complete/panUKBB_complete_clumping_250kb_r2_0.2.clumped | awk 'NF' > gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.250kb.0.2r2.prune.in
```
Map back from pos:id SNP IDs to RSIDs:
```bash
Rscript ${PU_SRC}/custom_rsid_map_post_clumpt.R ${GE}/panUKBB_complete_clumped_r2-0.2.250kb.0.2r2.prune.in ${GE}/panUKBB_complete_clumped_r2-0.2.pruned_rsids.txt
```

**12) Extract the matrix of summary statistics needed for GLEANR**

This requires as input the missingness_report.tsv and panUKBB_complete_clumped_r2-0.2.pruned_rsids.txt
```bash
snakemake --snakefile rules/extract_factorize.smk -j 1 gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.beta.tsv --configfile config/panUKBB_complete_config.yaml 
```

**13) Run GLEANR model selection**
```bash
	INTER_DIR=results/panUKBB_complete_41K_sklearn_eBIC
	mkdir -p ${INTER_DIR}
	#Because the data is large,I needed to run grid search on separate jobs (instead of running in one multi-threaded instance, which can be done natively by the gleanr packge)
	bash  run_scripts/GLEANER_panUKBB_41K_snps.build_sklearn-eBIC.sh #Identify the first 4 Kinit to test and build slurm scripts for that
```
Run model fitting at the first 4 points
```bash
sbatch final_gleaner_attempt_K102.slurm
sbatch final_gleaner_attempt_K136.slurm
sbatch final_gleaner_attempt_K34.slurm
sbatch final_gleaner_attempt_K68.slurm
```
Get the BIC for each of the above runs and nominate the next best K points to test
```bash
bash  run_scripts/GLEANER_panUKBB_41K_snps.score_sklearn-eBIC.sh
sbatch final_gleaner_attempt_K119.slurm
sbatch final_gleaner_attempt_K127.slurm
```
Get the BIC for all of the above runs and nominate the next best K points to test
```bash
bash  run_scripts/GLEANER_panUKBB_41K_snps.score_sklearn-eBIC.sh
sbatch final_gleaner_attempt_K134.slurm
sbatch final_gleaner_attempt_K130.slurm
```
Score all of the runs, and select the run with the lowest BIC
```bash
	bash  run_scripts/GLEANER_panUKBB_41K_snps.score_sklearn-eBIC.sh
```
 This nominates Kinit=134

**14) Run GLEANR model fitting step**
```bash
	sbatch final_gleaner_COMPLETION_RUN_K134.slurm
```
Move files into another directory for clean downstream analysis (optional)
```bash
RUN_DIR="results/panUKBB_complete_41K_sklearn_eBIC/"
TARGET_DIR="results/panUKBB_complete_41K_final/"
mkdir -p ${TARGET_DIR}
mkdir -p ${BASE_DIR}/panUKBB_complete_41K_final
mkdir -p results/panUKBB_complete_41K_final
cp ${RUN_DIR}/final_gleaner_attempt_K134_COMPLETIONlatent.factors.txt ${TARGET_DIR}/latent.factors.txt
cp ${RUN_DIR}/final_gleaner_attempt_K134_COMPLETIONlatent.loadings.txt ${TARGET_DIR}/latent.loadings.txt
cp gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.n.tsv gwas_extracts/panUKBB_complete_41K_final/panUKBB_complete_41K_final.n.tsv
cp gwas_extracts/panUKBB_complete/panUKBB_complete_clumped_r2-0.2.se.tsv gwas_extracts/panUKBB_complete_41K_final/panUKBB_complete_41K_final.se.tsv
cp ${RUN_DIR}/panUKBB_complete_41K_sklearn_eBIC/final_gleaner_attempt_K134_COMPLETION_final_dat.RData ${TARGET_DIR}/panUKBB_complete_41K_final_final_dat.RData
```
*In this run of gleanr there was a bug in PVE estimation resulting in a mixed order of factors. This has since been fixed, and all reported results are given in terms of the correctly ordered factors*

**15) Convert factors into Z-scores for S-LDSR analysis and run S-LDSC to get tissue-specific enrichments**
```bash
ml snakemake
rules/project_assess.smk
snakemake --snakefile rules/tissue_enrichment.smk   -j 5 results/panUKBB_complete_41K_final/NONE_ldsc_enrichment_Multi_tissue_chromatin/factor_global_fdr.heatmap.png --configfile config/panUKBB_complete_config.yaml --profile  profiles/rockfish/ 
```
**16) Perform per-factor gene set enrichment analysis**
```bash
snakemake --snakefile rules/interpret_factors.smk  results/panUKBB_complete_41K_final/top_elements_by_factor/top_fe/gene_factor_enrichments_unique/all_factors_evaluated.txt -j10
```
**17) Calculate statistics of factor genetic architecture**
```bash
snakemake --snakefile rules/interpret_factors.smk  results/panUKBB_complete_41K_final/selective_pressure/allFigs.RData -j1
```

**18) Generate figures and analysis from these results**

Key analyses shown in the manuscript can be found in the `manuscript_analyses/figure_scripts` directory:
- Traits selected and heatmap visualization of factors (Fig 3): `fig3_gleaner_heatmap.R`
- Evaluation of platelet related factors (Fig. 4): `fig_selectedion_blood_enrichment.R`
- Evaluation of factor architecture measures: `fig_4_factor_architecture.R`
